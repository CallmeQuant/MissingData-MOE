// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

// [[Rcpp::export]]
List critClustWrapper(NumericMatrix data, IntegerVector numExp, int k, std::string framework, std::string model, std::string crit="BIC", bool DA=false) {
    try {
        Environment base("package:base");
        Function dataframe = base["data.frame"];
        Function parseR = base["parse"];
        Function evalR = base["eval"];
        
        // Subset data columns based on numExp
        NumericMatrix dataAux(data.nrow(), numExp.size());
        for (int j = 0; j < numExp.size(); ++j) {
            dataAux(_, j) = data(_, numExp[j] - 1);
        }
        DataFrame df = dataframe(dataAux);

        if (framework == "Mclust") {
            Environment mclustEnv("package:mclust");
            Function MclustFunc = mclustEnv["Mclust"];
            Function hcFunction = mclustEnv["hc"];

            // List hcResult = hcFunction(Named("data") = df, Named("modelName") = "VVV");
            // List hcResult = hcFunction(Named("data") = df);
            // List mclustResult = MclustFunc(df, Named("G") = k,
            //                                Named("modelNames") = model,
            //                                Named("initialization") = List::create(Named("hcPairs") = hcResult));
            List mclustResult = MclustFunc(df, Named("G") = k,
                                            Named("modelNames") = model);
            return List::create(
                Named("criterionValue") = as<double>(mclustResult["bic"]),
                Named("criterion") = "BIC",
                Named("nbcluster") = as<int>(mclustResult["G"]),
                Named("model") = as<std::string>(mclustResult["modelName"]),
                Named("parameters") = mclustResult["parameters"],
                Named("proba") = mclustResult["z"],
                Named("partition") = mclustResult["classification"],
                Named("error") = "No error"
            );
        }
        else if (framework == "Rmixmod") {
            Environment Rmixmod("package:Rmixmod");
            Function RmixmodCluster = Rmixmod["mixmodCluster"];
            SEXP modelObject = evalR(parseR(Named("text") = model));
            
            S4 xem = RmixmodCluster(Named("data") = df,
                                    Named("nbCluster") = k,
                                    Named("models") = modelObject,
                                    Named("criterion") = crit);

            S4 bestResult = xem.slot("bestResult");

            return List::create(
                Named("criterionValue") = -as<double>(bestResult.slot("criterionValue")),
                Named("criterion") = as<std::string>(bestResult.slot("criterion")),
                Named("nbcluster") = as<int>(bestResult.slot("nbCluster")),
                Named("model") = as<std::string>(bestResult.slot("model")),
                Named("parameters") = bestResult.slot("parameters"),
                Named("proba") = bestResult.slot("proba"),
                Named("partition") = bestResult.slot("partition"),
                Named("error") = "No error"
            );
        }
        else { // MixAll
            Environment MixAll("package:MixAll");
            Function clusterDiagGaussian = MixAll["clusterDiagGaussian"];
            
            S4 xem = clusterDiagGaussian(Named("data") = df,
                                         Named("nbCluster") = k,
                                         Named("models") = model,
                                         Named("criterion") = crit,
                                         Named("nbCore") = 1);

            return List::create(
                Named("criterionValue") = -as<double>(xem.slot("criterion")),
                Named("criterion") = as<std::string>(xem.slot("criterionName")),
                Named("nbcluster") = as<int>(xem.slot("nbCluster")),
                Named("model") = model,
                Named("parameters") = xem.slot("component"),
                Named("proba") = xem.slot("tik"),
                Named("partition") = xem.slot("zi"),
                Named("error") = "No error"
            );
        }
    }
    catch (std::exception &ex) {
        Rcout << "Error: " << ex.what() << "\n";
        return List::create(Named("error") = ex.what());
    }
    catch (...) {
        Rcout << "Unknown error occurred.\n";
        return List::create(Named("error") = "Unknown error occurred.");
    }
}
