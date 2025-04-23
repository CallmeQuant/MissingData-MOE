// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <Rcpp.h>

// [[Rcpp::export]]
Rcpp::List mclustWrapper(Rcpp::NumericMatrix data, int k, std::string model) {
    try {
        // Load Mclust package and function
        Rcpp::Environment mclustEnv("package:mclust");
        Rcpp::Function MclustFunc = mclustEnv["Mclust"];
        
        // Call Mclust function
        Rcpp::CharacterVector modelVec = Rcpp::CharacterVector::create(model);
        Rcpp::List mclustResult = MclustFunc(Rcpp::_["data"] = data, 
                                            Rcpp::_["G"] = k, 
                                            Rcpp::_["modelNames"] = modelVec);
        
        // Extract values from Mclust result
        double bicVal = Rcpp::as<double>(mclustResult["bic"]);
        int nbCluster = Rcpp::as<int>(mclustResult["G"]);
        std::string modelName = Rcpp::as<std::string>(mclustResult["modelName"]);
        Rcpp::List parameters = mclustResult["parameters"];
        Rcpp::NumericMatrix proba = mclustResult["z"];
        Rcpp::IntegerVector partition = mclustResult["classification"];
        
        // Create an empty missing values data frame
        Rcpp::DataFrame missingVals = Rcpp::DataFrame::create(
            Rcpp::Named("row") = Rcpp::NumericVector(),
            Rcpp::Named("col") = Rcpp::NumericVector(),
            Rcpp::Named("value") = Rcpp::NumericVector()
        );
        
        // Construct result list
        Rcpp::List result = Rcpp::List::create(
            Rcpp::Named("criterionValue") = -bicVal,
            Rcpp::Named("criterion") = "BIC",
            Rcpp::Named("nbcluster") = nbCluster,
            Rcpp::Named("model") = modelName,
            Rcpp::Named("parameters") = parameters,
            Rcpp::Named("proba") = proba,
            Rcpp::Named("partition") = partition,
            Rcpp::Named("error") = "No error",
            Rcpp::Named("missingValues") = missingVals
        );
        return result;
    } catch (std::exception &ex) {
        Rcpp::Rcout << "Error in mclustWrapper: " << ex.what() << "\n";
        return Rcpp::List::create(Rcpp::Named("error") = ex.what());
    } catch (...) {
        Rcpp::Rcout << "Unknown error occurred in mclustWrapper.\n";
        return Rcpp::List::create(Rcpp::Named("error") = "Unknown error occurred in mclustWrapper");
    }
}
