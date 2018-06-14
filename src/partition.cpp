#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;
using namespace arma;

LogicalVector is_element(CharacterVector x, CharacterVector elements) {
  LogicalVector any_match(x.size()), element_match(elements.size());

  for(int i = 0; i < x.size(); ++i) {
   for(int k = 0; k < elements.size(); ++k) {
     element_match[k] = x[i] == elements[k];
   }
   any_match[i] = is_true(any(element_match));
  }
  return any_match;
}

IntegerVector which(LogicalVector x) {
  IntegerVector out;
  
  for(int i = 0; i < x.size(); i++) {
    if(x[i] == TRUE) {
      out.push_back(i);
     }
    }
  return out;
}

// [[Rcpp::export]]
NumericMatrix bind_empty_row(NumericMatrix x) {
  NumericVector row_vec(x.ncol());
  row_vec = rep(NumericVector::get_na(), x.ncol());
  NumericMatrix out(x.rows() + 1, x.cols());
  out = transpose(cbind(transpose(x), row_vec));
  return out;
}

// [[Rcpp::export]]
bool has_rownames(NumericMatrix x) {
  return !Rf_isNull(rownames(x));
}

// [[Rcpp::export]]
bool has_colnames(NumericMatrix x) {
  return !Rf_isNull(colnames(x));
}

// [[Rcpp::export]]
NumericMatrix subset_c(NumericMatrix x, IntegerVector subset) {
  NumericMatrix out(x.nrow(), subset.size());

  for (int i = 0; i < subset.size(); i++) {
    out.column(i) = x(_, subset[i]);
  }
  
  if (has_colnames(x)) {
    CharacterVector column_names, reduced_names;
    column_names = colnames(x);
    colnames(out) = column_names[subset];
  }
  
  if (has_rownames(x)) {
    CharacterVector row_names;
    row_names = rownames(x);
    rownames(out) = row_names;
  }
  
  return out;
}

// [[Rcpp::export]]
NumericMatrix drop(NumericMatrix x, IntegerVector subset) {

  NumericMatrix out(x.nrow(), x.ncol() - subset.size());
  LogicalVector should_drop(subset.size());
  
  int j = 0;
  for (int i = 0; i < x.ncol(); i++) {
    for (int k = 0; k < subset.size(); k++) {
      should_drop[k] = subset[k] == i;
    }
    
    if (is_false(any(should_drop))) {
      out.column(j) = x(_, i);
      j++;
    }
  }
  
  if (has_colnames(x)) {
    CharacterVector column_names;
    column_names = colnames(x);
    
    for (int k = 0; k < subset.size(); k++) {
      column_names.erase(subset[k]);
      subset = subset - 1;
    }
    
    colnames(out) = column_names;
  }
  
  if (has_rownames(x)) {
    CharacterVector row_names;
    row_names = rownames(x);
    rownames(out) = row_names;
  }
  
  return out;
}

// [[Rcpp::export]]
NumericMatrix drop_dist(NumericMatrix reduced_dist, IntegerVector distance_index) {
 IntegerVector distance_index_rows;
 distance_index_rows = clone(distance_index);
 reduced_dist = drop(reduced_dist, distance_index);
 reduced_dist = transpose(drop(transpose(reduced_dist), distance_index_rows));
 return reduced_dist;
}

// [[Rcpp::export]]
NumericVector icc_c(NumericMatrix x) {
  NumericVector icc, ms1, ms2, variance;
  int ncols = x.ncol(), nrows = x.nrow();
  NumericVector row_vec(ncols), rowmeans(nrows), long_means(nrows * ncols), 
    among, within;
  
  double matrix_mean = mean(x);
  
  for (int i = 0; i < nrows; ++i) {
    row_vec = x(i, _);
    rowmeans[i] = mean(row_vec);
    row_vec = row_vec - rowmeans[i];
    for (NumericVector::iterator it = row_vec.begin(); it != row_vec.end(); ++it) {
      within.insert(within.end(), *it);
    }
  }
  
  long_means = rep(rowmeans, ncols);
  among = pow(long_means - matrix_mean, 2);
  within = pow(within, 2);
  ms1 = sum(among) / (nrows - 1); 
  ms2 = sum(within) / (nrows * (ncols - 1));
  variance = (ms1 - ms2) / ncols;
  icc = variance / (variance + ms2);
  
  return icc;
}

// [[Rcpp::export]]
NumericVector scale_rowmeans(NumericMatrix x) {
  
  NumericVector row_vec(x.ncol()), out(x.nrow());
  
  for(int i = 0; i < x.nrow(); ++i) {
    row_vec = x(i, _);
    out[i] = mean(noNA(row_vec));
  }
  
  return (out - mean(out)) / sd(out);
}

// [[Rcpp::export]]
List ICC_c(NumericMatrix x) {
  return List::create(Named("ICC", icc_c(x)),
                      Named("row_means", scale_rowmeans(x)));
}


class SortRanks {
private:
    const Rcpp::NumericVector& ref;

    bool is_na(double x) const 
    {
        return Rcpp::traits::is_na<REALSXP>(x);    
    }

public:
    SortRanks(const Rcpp::NumericVector& ref_)
        : ref(ref_)
    {}

    bool operator()(const int ilhs, const int irhs) const
    {
        double lhs = ref[ilhs], rhs = ref[irhs]; 
        if (is_na(lhs)) return false;
        if (is_na(rhs)) return true;
        return lhs < rhs;
    }
};

// [[Rcpp::export]]
Rcpp::NumericVector rank_c(Rcpp::NumericVector x)
{
    int vec_size = x.size();
    Rcpp::IntegerVector ranks = seq(0, vec_size - 1);
    std::sort(ranks.begin(), ranks.end(), SortRanks(x));

    Rcpp::NumericVector avg_ranks(vec_size);
    for (int n, i = 0; i < vec_size; i += n) {
        n = 1;
        while (i + n < vec_size && x[ranks[i]] == x[ranks[i + n]]) ++n;
        for (int k = 0; k < n; k++) {
            avg_ranks[ranks[i + k]] = i + (n + 1) / 2.0;
        }
    }

    return avg_ranks;
}

// [[Rcpp::export]]
arma::mat apply_rank(arma::mat x) {
  NumericVector temp_vec;
  NumericMatrix out(x.n_rows, x.n_cols);
  NumericMatrix y = wrap(x);
  
  for(int i = 0; i < y.ncol(); ++i) {
    temp_vec = y(_, i);
    out(_, i) = rank_c(temp_vec);
  }
  
  return as<arma::mat>(out);
}

// [[Rcpp::export]]
arma::mat corr_c_mat(arma::mat x) {
  return arma::cor(x);
}

// [[Rcpp::export]]
arma::mat corr_c_2mat(arma::mat x, arma::mat y) {
  return arma::cor(x, y);
}

// [[Rcpp::export]]
double corr_c_2vec(arma::vec x, arma::vec y) {
  arma::mat out;
  out = arma::cor(x, y);
  return out(0);
}

// [[Rcpp::export]]
NumericMatrix pearson_distance(NumericMatrix x, NumericMatrix y) {
  return wrap(1 - corr_c_2mat(as<arma::mat>(x), as<arma::mat>(y)));
}

// [[Rcpp::export]]
NumericMatrix spearman_distance(NumericMatrix x, NumericMatrix y) {
  return wrap(1 - corr_c_2mat(apply_rank(as<arma::mat>(x)), apply_rank(as<arma::mat>(y))));
}

// [[Rcpp::export]]
List pca_c(arma::mat x) {
  NumericMatrix rcpp_x = wrap(x);
  NumericVector col_vec;
  
  // scale each variable
  for(int i = 0; i < rcpp_x.ncol(); ++i) {
    col_vec = rcpp_x(_, i);
    rcpp_x(_, i) = (col_vec - mean(noNA(col_vec))) / sd(noNA(col_vec));
  }
  
 x = as<arma::mat>(rcpp_x);
  
 mat coeff;
 mat score;
 vec latent, pc1;
 double pct_var;

 arma::princomp(coeff, score, latent, x);
 
 pc1 = score.col(0);
 pc1 = (pc1 - mean(pc1)) / stddev(pc1);
 pct_var = latent[0] / sum(latent);
 
 return List::create(pct_var, pc1);
}

// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;

// [[Rcpp::export]]
List minR2_c(NumericMatrix x) {
  
  NumericVector row_vec(x.ncol()), row_means(x.nrow());
  NumericMatrix cors(x.ncol() + 1), x_means(x.ncol() + 1);
  double r2;
  
  for(int i = 0; i < x.nrow(); ++i) {
    row_vec = x(i, _);
    row_means[i] = mean(noNA(row_vec));
  }
  
  row_means = (row_means - mean(row_means)) / sd(row_means);
  x_means = cbind(row_means, x);
  cors = wrap(arma::cor(as<arma::mat>(x_means)));
  
	for(int i = 0; i < cors.nrow(); ++i) cors(i, _) = pow(cors(i, _), 2);
  
  r2 = min(cors(_, 1));
  
  return List::create(r2, row_means);
}

// [[Rcpp::export]]
NumericMatrix update_dist(NumericMatrix reduced_dist, CharacterVector cluster_nm,  
                          CharacterVector clust_var_nms, NumericMatrix reduced_data, 
                          std::string dist_type) {
  IntegerVector cluster_index, distance_index;
  NumericMatrix cluster_var, other_vars, cluster_dist;
  cluster_index = which(is_element(colnames(reduced_data), cluster_nm));
                           
	cluster_var = subset_c(reduced_data, cluster_index);
	other_vars = drop(reduced_data, cluster_index);

	//nms = colnames(dat.r)[colnames(dat.r) %nin% cluster.nm]
	//dat.r.tmp = as.data.frame( dat.r[, nms ] )
	//names( dat.r.tmp ) = nms
	// compute distances between myvar and other variables
	//if( !is.na(dim(dat.r.tmp)[1]) ){
		if (dist_type == "s") {
			cluster_dist = spearman_distance(other_vars, cluster_var);
		} else {
			cluster_dist = pearson_distance(other_vars, cluster_var);
		}
		//colnames(cluster_dist) = colnames(other_vars);
	//}

	// remove cluster component variables from distance matrix
  distance_index = which(is_element(colnames(reduced_dist), clust_var_nms));
	reduced_dist = drop_dist(reduced_dist, distance_index);

	//if( !is.null(dim(dist.r)[1])){
	  CharacterVector reduced_colnames;
	  reduced_colnames = colnames(reduced_dist);
	  reduced_colnames.push_back(cluster_nm[0]);
	 
	  NumericMatrix out_dist, bound_dist;
    bound_dist = cbind(reduced_dist, cluster_dist);
    

//     CharacterVector reduced_rownames;
// 		reduced_rownames = rownames(reduced_dist);
// 	  reduced_rownames.push_back(cluster_nm[0]);
		out_dist = bind_empty_row(bound_dist);
		colnames(out_dist) = reduced_colnames;
		rownames(out_dist) = reduced_colnames;

	//}
	
	return out_dist;
} 


// // [[Rcpp::export]]
// List assign_clusters(IntegerVector index, NumericMatrix reduced_dist, 
//                      NumericMatrix reduced_data_r, NumericMatrix data_r, 
//                      double pct_var, DataFrame clusters_r, int cluster_ind, 
//                      std::string method, std::string dist_type, std::string new_var){
// 	bool success = false;
//   NumericMatrix data, reduced_data;
//   DataFrame clusters = clone(clusters_r);
//   data = clone(data_r);
//   reduced_data = clone(reduced_data_r);
// 	CharacterVector clust_colnames, original_vars, row_names, clust_var_nms, clust_var_nms_raw;
//   IntegerVector var_ind;
//   clust_colnames = colnames(reduced_data);
//   original_vars = clusters.attr("row.names");
//   row_names = clone(original_vars);
//   clust_var_nms = clust_colnames[index - 1];
//   
// 	clust_var_nms_raw = original_vars[is_element(clusters["cluster"], clust_var_nms)];
// 	std::string cluster_nm = new_var.append(std::to_string(cluster_ind));
// 	
// 	// 2. compute summary and test against pct.var
// 	
// 	var_ind = which(is_element(original_vars, clust_var_nms_raw));
// 	List tmp_svec; 
// 	NumericMatrix x = subset_c(data, var_ind);
// 	
// 	if (method == "ICC") {
// 	  tmp_svec = ICC_c(x);
// 	} else if (method == "MI") {
// 	  Function MI("MI");
// 	  tmp_svec = MI(x);
// 	} else if (method == "minR2") {
// 	  tmp_svec = minR2_c(x);
// 	} else {
// 	  tmp_svec = pca_c(as<arma::mat>(x));
// 	}
//                     
// 
// 	// 3. update, dist.r, dat.r, clusters, cluster.new.names
// 	double clstr_pct_var;
// 	NumericVector clstr_data;
// 	clstr_pct_var = tmp_svec[0];
// 	clstr_data = tmp_svec[1];
// 	
//   CharacterVector cluster_cluster;
//   NumericVector cluster_pct_var;
//   LogicalVector cluster_index;
//   DataFrame mappings;
//   NumericMatrix out_data, out_dist;
//   
// 	if( clstr_pct_var >= pct_var ){
// 			success = true;
// 
// 	    // update cluster mappings
// 
// 	    cluster_cluster = clusters["cluster"];
// 	    cluster_pct_var = clusters["pct.var"];
// 	    cluster_index = is_element(cluster_cluster, clust_var_nms_raw);
// 
// 	    for (int i = 0; i < cluster_index.size(); i++) {
// 	      if (cluster_index[i]) {
// 	        cluster_cluster[i] = cluster_nm;
// 	        cluster_pct_var[i] = clstr_pct_var;
//   	      }
//   	    }
// 
// 			mappings = DataFrame::create(_["cluster"]= cluster_cluster,
//                      _["pct.var"] = cluster_pct_var);
//       mappings.attr("row.names") = row_names;
// 			// increment to form unique name for next formed cluster
// 			cluster_ind++;
// 
// 			// add to reduced data and update names
// 			out_data = cbind(reduced_data, clstr_data); // add summary variable to dat.r
// 			CharacterVector new_colnames;
// 			new_colnames = colnames(reduced_data);
// 			new_colnames.push_back(cluster_nm);
// 			colnames(out_data) = new_colnames;
// 			out_data = drop(out_data, index);
// 			//  aa = !is.null(dim(dat.r)[1])
// 			//if( aa )
// 			
// 		out_dist = update_dist(reduced_dist, cluster_nm, clust_var_nms, out_data, dist_type);
// 	} else {
// 	  // set distance to NA to avoid testing again
// 		reduced_dist(index[0], index[1]) = NumericMatrix::get_na();
// 	  out_data = reduced_data;
// 	  out_dist = reduced_dist;
// 	}
// 	
// 	return List::create(mappings, cluster_ind, out_dist, out_data, success);
// }


// [[Rcpp::export]]
List assign_clusters(IntegerVector index_r, NumericMatrix reduced_dist_r, 
                     NumericMatrix reduced_data_r, NumericMatrix data_r, 
                     double pct_var, DataFrame clusters_r, int cluster_ind, 
                     std::string method, std::string dist_type, std::string new_var){
	bool success = false;
  NumericMatrix data, reduced_data, reduced_dist;
  DataFrame clusters = clone(clusters_r);
  data = clone(data_r);
  reduced_data = clone(reduced_data_r);
  reduced_dist = clone(reduced_dist_r);
	CharacterVector clust_colnames, original_vars, row_names, clust_var_nms, clust_var_nms_raw;
  IntegerVector var_ind, index;
  index = clone(index_r);
  clust_colnames = colnames(reduced_data);
  original_vars = clusters.attr("row.names");
  row_names = clone(original_vars);
  clust_var_nms = clust_colnames[index - 1];
  
	clust_var_nms_raw = original_vars[is_element(row_names, clust_var_nms)];
	std::string cluster_nm = new_var.append(std::to_string(cluster_ind));
	
	// 2. compute summary and test against pct.var
	
	var_ind = which(is_element(original_vars, clust_var_nms_raw));
	List tmp_svec; 
	NumericMatrix x = subset_c(data, var_ind);
	
	if (method == "ICC") {
	  tmp_svec = ICC_c(x);
	} else if (method == "MI") {
	  Environment pkg = Environment::namespace_env("partition");
    Function MI = pkg["MI"];
	  tmp_svec = MI(x);
	} else if (method == "minR2") {
	  tmp_svec = minR2_c(x);
	} else {
	  tmp_svec = pca_c(as<arma::mat>(x));
	}
                    

	// 3. update, dist.r, dat.r, clusters, cluster.new.names
	double clstr_pct_var;
	NumericVector clstr_data;
	clstr_pct_var = tmp_svec[0];
	clstr_data = tmp_svec[1];
	
  CharacterVector cluster_cluster;
  NumericVector cluster_pct_var;
  LogicalVector cluster_index;
  DataFrame mappings;
  NumericMatrix out_data, out_dist;
  
	if( clstr_pct_var >= pct_var ){
			success = true;

	    // update cluster mappings

	    cluster_cluster = clusters["cluster"];
	    cluster_pct_var = clusters["pct.var"];
	    cluster_index = is_element(row_names, clust_var_nms_raw);

	    for (int i = 0; i < cluster_index.size(); i++) {
	      if (cluster_index[i]) {
	        cluster_cluster[i] = cluster_nm;
	        cluster_pct_var[i] = clstr_pct_var;
  	      }
  	    }

			mappings = DataFrame::create(_["cluster"]= cluster_cluster,
                     _["pct.var"] = cluster_pct_var);
      mappings.attr("row.names") = row_names;
			// increment to form unique name for next formed cluster
			cluster_ind++;

			// add to reduced data and update names
			// CharacterVector row_incides;
			// row_incides = rownames(reduced_data);
			out_data = cbind(reduced_data, clstr_data); // add summary variable to dat.r
			CharacterVector new_colnames;
			new_colnames = colnames(reduced_data);
			new_colnames.push_back(cluster_nm);
			colnames(out_data) = new_colnames;
			out_data = drop(out_data, index - 1);
			
			//rownames(out_data) = row_incides;
			//  aa = !is.null(dim(dat.r)[1])
			//if( aa )
			
		out_dist = update_dist(reduced_dist, cluster_nm, clust_var_nms, out_data, dist_type);
	} else {
	  // set distance to NA to avoid testing again
		reduced_dist(index[0] - 1, index[1] - 1) = NumericMatrix::get_na();
	  out_data = reduced_data;
	  mappings = clusters;
	  out_dist = reduced_dist;
	}
	
	return List::create(mappings, cluster_ind, out_dist, out_data, success);
}