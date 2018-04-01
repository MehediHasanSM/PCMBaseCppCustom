/*
 *  Rcpp.cpp
 *  PCMBaseCpp
 *
 * Copyright 2017 Venelin Mitov
 *
 * This file is part of PCMBaseCpp: A C++ backend for calculating the likelihood of phylogenetic comparative models.
 *
 * PCMBaseCpp is free software: you can redistribute it and/or modify
 * it under the terms of version 3 of the GNU General Public License as
 * published by the Free Software Foundation.
 *
 * PCMBaseCpp is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public
 * License along with PCMBaseCpp.  If not, see
 * <http://www.gnu.org/licenses/>.
 *
 * @author Venelin Mitov
 */
#include <RcppArmadillo.h>
#include <R_ext/Rdynload.h>

#include<vector>
#include<string>
#include<unordered_map>
#include<sstream>
#include<iostream>

#include "QuadraticPolynomialWhite.h"
#include "QuadraticPolynomialBM.h"
#include "QuadraticPolynomialOU.h"
#include "QuadraticPolynomialJOU.h"
#include "QuadraticPolynomialTwoSpeedOU.h"
#include "QuadraticPolynomialMRG.h"

// [[Rcpp::plugins("cpp11")]]
// [[Rcpp::plugins(openmp)]]
// [[Rcpp::depends(RcppArmadillo)]]

// BEGIN: Needed for r-devel (R 3.4)
void R_init_PCMBaseCpp(DllInfo *info) {
  /* Register routines, allocate resources. */
  R_registerRoutines(info, NULL, NULL, NULL, NULL);
  R_useDynamicSymbols(info, TRUE);
}

void R_unload_PCMBaseCpp(DllInfo *info) {
  /* Release resources. */
}
// END Needed for r-devel (R 3.4)



// [[Rcpp::plugins("cpp11")]]
// [[Rcpp::plugins(openmp)]]
// [[Rcpp::depends(RcppArmadillo)]]

using namespace PCMBaseCpp;
using namespace std;

QuadraticPolynomialWhite* CreateQuadraticPolynomialWhite(
    arma::mat const&X, 
    Rcpp::List const& tree, 
    Rcpp::List const& model,
    Rcpp::List const& metaInfo) { 
  
  double threshold_SV = static_cast<double>(metaInfo["PCMBase.Threshold.SV"]);
  double threshold_skip_singular = static_cast<double>(metaInfo["PCMBase.Threshold.Skip.Singular"]);
  double threshold_Lambda_ij = static_cast<double>(metaInfo["PCMBase.Threshold.Lambda_ij"]);
  bool internal_pc_full = static_cast<int>(metaInfo["PCMBase.Internal.PC.Full"]);
  bool skip_singular = static_cast<int>(metaInfo["PCMBase.Skip.Singular"]);
  
  
  arma::imat Pc(X.n_rows, X.n_cols, arma::fill::ones);
  
  for(arma::uword i = 0; i < X.n_rows; ++i)
    for(arma::uword j = 0; j < X.n_cols; ++j) {
      Pc(i,j) = static_cast<arma::imat::value_type>(arma::is_finite(X(i,j)));
    }
    
    arma::umat branches = tree["edge"];
  splittree::uvec br_0 = arma::conv_to<splittree::uvec>::from(branches.col(0));
  splittree::uvec br_1 = arma::conv_to<splittree::uvec>::from(branches.col(1));
  splittree::vec t = Rcpp::as<splittree::vec>(tree["edge.length"]);
  
  using namespace std;
  
  uint RModel = Rcpp::as<uint>(metaInfo["RModel"]);
  
  vector<arma::uword> regimes = Rcpp::as<vector<arma::uword>>(metaInfo["r"]);
  
  if(regimes.size() != branches.n_rows) {
    ostringstream os;
    os<<"ERR:03801:PCMBaseCpp:Rcpp.cpp:CreateQuadraticPolynomialWhite:: The slot r in metaInfo has different length ("<<regimes.size()<<
      ") than the number of edges ("<<branches.n_rows<<").";
    throw logic_error(os.str());
  }
  
  splittree::uint num_tips = Rcpp::as<Rcpp::CharacterVector>(tree["tip.label"]).size();
  splittree::uvec node_names = splittree::Seq(static_cast<splittree::uint>(1), num_tips);
  
  vector<typename QuadraticPolynomialWhite::LengthType> lengths(branches.n_rows);
  
  for(arma::uword i = 0; i < branches.n_rows; ++i) {
    lengths[i].length_ = t[i];
    lengths[i].regime_ = regimes[i] - 1;
  }
  
  
  if(threshold_SV <= 0) {
    ostringstream os;
    os<<"ERR:03802:PCMBaseCpp:Rcpp.cpp:CreateQuadraticPolynomialWhite:: The argument threshold_SV should be positive real number.";
    throw invalid_argument(os.str());
  }
  
  typename QuadraticPolynomialWhite::DataType data(
      node_names, X, Pc, internal_pc_full, RModel, std::vector<std::string>(), 
      threshold_SV, threshold_skip_singular, skip_singular,
      threshold_Lambda_ij);
  
  return new QuadraticPolynomialWhite(br_0, br_1, lengths, data);
}

  RCPP_EXPOSED_CLASS_NODECL(QuadraticPolynomialWhite::TreeType)
  RCPP_EXPOSED_CLASS_NODECL(QuadraticPolynomialWhite::TraversalSpecificationType)
  RCPP_EXPOSED_CLASS_NODECL(QuadraticPolynomialWhite::AlgorithmType)
  
  RCPP_MODULE(QuadraticPolynomialWhite) {
    Rcpp::class_<QuadraticPolynomialWhite::TreeType::Tree> ( "QuadraticPolynomialWhite_Tree" )
    .property("num_nodes", &QuadraticPolynomialWhite::TreeType::Tree::num_nodes )
    .property("num_tips", &QuadraticPolynomialWhite::TreeType::Tree::num_tips )
    .method("FindNodeWithId", &QuadraticPolynomialWhite::TreeType::Tree::FindNodeWithId )
    .method("FindIdOfNode", &QuadraticPolynomialWhite::TreeType::Tree::FindIdOfNode )
    .method("FindIdOfParent", &QuadraticPolynomialWhite::TreeType::Tree::FindIdOfParent )
    .method("OrderNodes", &QuadraticPolynomialWhite::TreeType::Tree::OrderNodes )
    ;
    Rcpp::class_<QuadraticPolynomialWhite::TreeType>( "QuadraticPolynomialWhite_OrderedTree" )
      .derives<QuadraticPolynomialWhite::TreeType::Tree> ( "QuadraticPolynomialWhite_Tree" )
      .method("RangeIdPruneNode", &QuadraticPolynomialWhite::TreeType::RangeIdPruneNode )
      .method("RangeIdVisitNode", &QuadraticPolynomialWhite::TreeType::RangeIdVisitNode )
      .property("num_levels", &QuadraticPolynomialWhite::TreeType::num_levels )
      .property("ranges_id_visit", &QuadraticPolynomialWhite::TreeType::ranges_id_visit )
      .property("ranges_id_prune", &QuadraticPolynomialWhite::TreeType::ranges_id_prune )
    ;
    Rcpp::class_<QuadraticPolynomialWhite::AlgorithmType::ParentType>( "QuadraticPolynomialWhite_TraversalAlgorithm" )
      .property( "VersionOPENMP", &QuadraticPolynomialWhite::AlgorithmType::ParentType::VersionOPENMP )
      .property( "num_threads", &QuadraticPolynomialWhite::AlgorithmType::num_threads )
    ;
    Rcpp::class_<QuadraticPolynomialWhite::AlgorithmType> ( "QuadraticPolynomialWhite_ParallelPruning" )
      .derives<QuadraticPolynomialWhite::AlgorithmType::ParentType>( "QuadraticPolynomialWhite_TraversalAlgorithm" )
      .method( "ModeAutoStep", &QuadraticPolynomialWhite::AlgorithmType::ModeAutoStep )
      .property( "ModeAutoCurrent", &QuadraticPolynomialWhite::AlgorithmType::ModeAutoCurrent )
      .property( "IsTuning", &QuadraticPolynomialWhite::AlgorithmType::IsTuning )
      .property( "min_size_chunk_visit", &QuadraticPolynomialWhite::AlgorithmType::min_size_chunk_visit )
      .property( "min_size_chunk_prune", &QuadraticPolynomialWhite::AlgorithmType::min_size_chunk_prune )
      .property( "durations_tuning", &QuadraticPolynomialWhite::AlgorithmType::durations_tuning )
      .property( "fastest_step_tuning", &QuadraticPolynomialWhite::AlgorithmType::fastest_step_tuning )
    ;
    Rcpp::class_<QuadraticPolynomialWhite::TraversalSpecificationType::BaseType> ( "QuadraticPolynomial_PrunignSpec" )
      .field( "A", &QuadraticPolynomialWhite::TraversalSpecificationType::BaseType::A)
      .field( "b", &QuadraticPolynomialWhite::TraversalSpecificationType::BaseType::b )
      .field( "C", &QuadraticPolynomialWhite::TraversalSpecificationType::BaseType::C )
      .field( "d", &QuadraticPolynomialWhite::TraversalSpecificationType::BaseType::d )
      .field( "E", &QuadraticPolynomialWhite::TraversalSpecificationType::BaseType::E )
      .field( "f", &QuadraticPolynomialWhite::TraversalSpecificationType::BaseType::f )
      .field( "L", &QuadraticPolynomialWhite::TraversalSpecificationType::BaseType::L )
      .field( "m", &QuadraticPolynomialWhite::TraversalSpecificationType::BaseType::m )
      .field( "r", &QuadraticPolynomialWhite::TraversalSpecificationType::BaseType::r )
      .field( "V", &QuadraticPolynomialWhite::TraversalSpecificationType::BaseType::V )
      .field( "V_1", &QuadraticPolynomialWhite::TraversalSpecificationType::BaseType::V_1 )
    ;
    Rcpp::class_<QuadraticPolynomialWhite::TraversalSpecificationType> ( "QuadraticPolynomialWhite_PruningSpec" )
      .derives<QuadraticPolynomialWhite::TraversalSpecificationType::BaseType>( "QuadraticPolynomial_PrunignSpec" )
    ;
    Rcpp::class_<QuadraticPolynomialWhite>( "QuadraticPolynomialWhite" )
      .factory<arma::mat const&, Rcpp::List const&, Rcpp::List const&>(&CreateQuadraticPolynomialWhite)
      .method( "TraverseTree", &QuadraticPolynomialWhite::TraverseTree )
      .property( "tree", &QuadraticPolynomialWhite::tree )
      .property( "spec", &QuadraticPolynomialWhite::spec )
      .property( "algorithm", &QuadraticPolynomialWhite::algorithm )
    ;
  }

QuadraticPolynomialBM* CreateQuadraticPolynomialBM(
    arma::mat const&X, 
    Rcpp::List const& tree, 
    Rcpp::List const& model,
    Rcpp::List const& metaInfo) { 
    
  double threshold_SV = static_cast<double>(metaInfo["PCMBase.Threshold.SV"]);
  double threshold_skip_singular = static_cast<double>(metaInfo["PCMBase.Threshold.Skip.Singular"]);
  double threshold_Lambda_ij = static_cast<double>(metaInfo["PCMBase.Threshold.Lambda_ij"]);
  bool internal_pc_full = static_cast<int>(metaInfo["PCMBase.Internal.PC.Full"]);
  bool skip_singular = static_cast<int>(metaInfo["PCMBase.Skip.Singular"]);
  
  
  arma::imat Pc(X.n_rows, X.n_cols, arma::fill::ones);
  
  for(arma::uword i = 0; i < X.n_rows; ++i)
    for(arma::uword j = 0; j < X.n_cols; ++j) {
      Pc(i,j) = static_cast<arma::imat::value_type>(arma::is_finite(X(i,j)));
    }
    
  arma::umat branches = tree["edge"];
  splittree::uvec br_0 = arma::conv_to<splittree::uvec>::from(branches.col(0));
  splittree::uvec br_1 = arma::conv_to<splittree::uvec>::from(branches.col(1));
  splittree::vec t = Rcpp::as<splittree::vec>(tree["edge.length"]);
  
  using namespace std;
  
  uint RModel = Rcpp::as<uint>(metaInfo["RModel"]);
  
  vector<arma::uword> regimes = Rcpp::as<vector<arma::uword>>(metaInfo["r"]);
  
  if(regimes.size() != branches.n_rows) {
    ostringstream os;
    os<<"ERR:03801:PCMBaseCpp:Rcpp.cpp:CreateQuadraticPolynomialBM:: The slot r in metaInfo has different length ("<<regimes.size()<<
      ") than the number of edges ("<<branches.n_rows<<").";
    throw logic_error(os.str());
  }
  
  splittree::uint num_tips = Rcpp::as<Rcpp::CharacterVector>(tree["tip.label"]).size();
  splittree::uvec node_names = splittree::Seq(static_cast<splittree::uint>(1), num_tips);
  
  vector<typename QuadraticPolynomialBM::LengthType> lengths(branches.n_rows);
  
  for(arma::uword i = 0; i < branches.n_rows; ++i) {
    lengths[i].length_ = t[i];
    lengths[i].regime_ = regimes[i] - 1;
  }
  
  
  if(threshold_SV <= 0) {
    ostringstream os;
    os<<"ERR:03802:PCMBaseCpp:Rcpp.cpp:CreateQuadraticPolynomialBM:: The argument threshold_SV should be positive real number.";
    throw invalid_argument(os.str());
  }
  
  typename QuadraticPolynomialBM::DataType data(
      node_names, X, Pc, internal_pc_full, RModel, std::vector<std::string>(), 
      threshold_SV, threshold_skip_singular, skip_singular,
      threshold_Lambda_ij);
  
  return new QuadraticPolynomialBM(br_0, br_1, lengths, data);
}

//RCPP_EXPOSED_CLASS_NODECL(QuadraticPolynomialBM::TreeType)
RCPP_EXPOSED_CLASS_NODECL(QuadraticPolynomialBM::TraversalSpecificationType)
RCPP_EXPOSED_CLASS_NODECL(QuadraticPolynomialBM::AlgorithmType)
  
RCPP_MODULE(QuadraticPolynomialBM) {
  Rcpp::class_<QuadraticPolynomialBM::TreeType::Tree> ( "QuadraticPolynomialBM_Tree" )
  .property("num_nodes", &QuadraticPolynomialBM::TreeType::Tree::num_nodes )
  .property("num_tips", &QuadraticPolynomialBM::TreeType::Tree::num_tips )
  .method("FindNodeWithId", &QuadraticPolynomialBM::TreeType::Tree::FindNodeWithId )
  .method("FindIdOfNode", &QuadraticPolynomialBM::TreeType::Tree::FindIdOfNode )
  .method("FindIdOfParent", &QuadraticPolynomialBM::TreeType::Tree::FindIdOfParent )
  .method("OrderNodes", &QuadraticPolynomialBM::TreeType::Tree::OrderNodes )
  ;
  Rcpp::class_<QuadraticPolynomialBM::TreeType>( "QuadraticPolynomialBM_OrderedTree" )
    .derives<QuadraticPolynomialBM::TreeType::Tree> ( "QuadraticPolynomialBM_Tree" )
    .method("RangeIdPruneNode", &QuadraticPolynomialBM::TreeType::RangeIdPruneNode )
    .method("RangeIdVisitNode", &QuadraticPolynomialBM::TreeType::RangeIdVisitNode )
    .property("num_levels", &QuadraticPolynomialBM::TreeType::num_levels )
    .property("ranges_id_visit", &QuadraticPolynomialBM::TreeType::ranges_id_visit )
    .property("ranges_id_prune", &QuadraticPolynomialBM::TreeType::ranges_id_prune )
  ;
  Rcpp::class_<QuadraticPolynomialBM::AlgorithmType::ParentType>( "QuadraticPolynomialBM_TraversalAlgorithm" )
    .property( "VersionOPENMP", &QuadraticPolynomialBM::AlgorithmType::ParentType::VersionOPENMP )
    .property( "num_threads", &QuadraticPolynomialBM::AlgorithmType::num_threads )
  ;
  Rcpp::class_<QuadraticPolynomialBM::AlgorithmType> ( "QuadraticPolynomialBM_ParallelPruning" )
    .derives<QuadraticPolynomialBM::AlgorithmType::ParentType>( "QuadraticPolynomialBM_TraversalAlgorithm" )
    .method( "ModeAutoStep", &QuadraticPolynomialBM::AlgorithmType::ModeAutoStep )
    .property( "ModeAutoCurrent", &QuadraticPolynomialBM::AlgorithmType::ModeAutoCurrent )
    .property( "IsTuning", &QuadraticPolynomialBM::AlgorithmType::IsTuning )
    .property( "min_size_chunk_visit", &QuadraticPolynomialBM::AlgorithmType::min_size_chunk_visit )
    .property( "min_size_chunk_prune", &QuadraticPolynomialBM::AlgorithmType::min_size_chunk_prune )
    .property( "durations_tuning", &QuadraticPolynomialBM::AlgorithmType::durations_tuning )
    .property( "fastest_step_tuning", &QuadraticPolynomialBM::AlgorithmType::fastest_step_tuning )
  ;
  Rcpp::class_<QuadraticPolynomialBM::TraversalSpecificationType::BaseType> ( "QuadraticPolynomial_PrunignSpec" )
    .field( "A", &QuadraticPolynomialBM::TraversalSpecificationType::BaseType::A)
    .field( "b", &QuadraticPolynomialBM::TraversalSpecificationType::BaseType::b )
    .field( "C", &QuadraticPolynomialBM::TraversalSpecificationType::BaseType::C )
    .field( "d", &QuadraticPolynomialBM::TraversalSpecificationType::BaseType::d )
    .field( "E", &QuadraticPolynomialBM::TraversalSpecificationType::BaseType::E )
    .field( "f", &QuadraticPolynomialBM::TraversalSpecificationType::BaseType::f )
    .field( "L", &QuadraticPolynomialBM::TraversalSpecificationType::BaseType::L )
    .field( "m", &QuadraticPolynomialBM::TraversalSpecificationType::BaseType::m )
    .field( "r", &QuadraticPolynomialBM::TraversalSpecificationType::BaseType::r )
    .field( "V", &QuadraticPolynomialBM::TraversalSpecificationType::BaseType::V )
    .field( "V_1", &QuadraticPolynomialBM::TraversalSpecificationType::BaseType::V_1 )
  ;
  Rcpp::class_<QuadraticPolynomialBM::TraversalSpecificationType> ( "QuadraticPolynomialBM_PruningSpec" )
    .derives<QuadraticPolynomialBM::TraversalSpecificationType::BaseType>( "QuadraticPolynomial_PrunignSpec" )
  ;
  Rcpp::class_<QuadraticPolynomialBM>( "QuadraticPolynomialBM" )
    .factory<arma::mat const&, Rcpp::List const&, Rcpp::List const&>(&CreateQuadraticPolynomialBM)
    .method( "TraverseTree", &QuadraticPolynomialBM::TraverseTree )
    .property( "tree", &QuadraticPolynomialBM::tree )
    .property( "spec", &QuadraticPolynomialBM::spec )
    .property( "algorithm", &QuadraticPolynomialBM::algorithm )
  ;
}

QuadraticPolynomialOU* CreateQuadraticPolynomialOU(
    arma::mat const&X, 
    Rcpp::List const& tree, 
    Rcpp::List const& model,
    Rcpp::List const& metaInfo) {
  
  double threshold_SV = static_cast<double>(metaInfo["PCMBase.Threshold.SV"]);
  double threshold_skip_singular = static_cast<double>(metaInfo["PCMBase.Threshold.Skip.Singular"]);
  double threshold_Lambda_ij = static_cast<double>(metaInfo["PCMBase.Threshold.Lambda_ij"]);
  bool internal_pc_full = static_cast<int>(metaInfo["PCMBase.Internal.PC.Full"]);
  bool skip_singular = static_cast<int>(metaInfo["PCMBase.Skip.Singular"]);
  
  arma::imat Pc(X.n_rows, X.n_cols, arma::fill::ones);
  
  for(arma::uword i = 0; i < X.n_rows; ++i)
    for(arma::uword j = 0; j < X.n_cols; ++j) {
      Pc(i,j) = static_cast<arma::imat::value_type>(arma::is_finite(X(i,j)));
    }
    
  arma::umat branches = tree["edge"];
  splittree::uvec br_0 = arma::conv_to<splittree::uvec>::from(branches.col(0));
  splittree::uvec br_1 = arma::conv_to<splittree::uvec>::from(branches.col(1));
  splittree::vec t = Rcpp::as<splittree::vec>(tree["edge.length"]);
  
  using namespace std;

  uint RModel = Rcpp::as<uint>(metaInfo["RModel"]);
  
  vector<arma::uword> regimes = Rcpp::as<vector<arma::uword>>(metaInfo["r"]);
  

  if(regimes.size() != branches.n_rows) {
    ostringstream os;
    os<<"ERR:03811:PCMBaseCpp:Rcpp.cpp:CreateQuadraticPolynomialOU:: The slot r in metaInfo has different length ("<<regimes.size()<<
      ") than the number of edges ("<<branches.n_rows<<").";
    throw logic_error(os.str());
  }

  splittree::uint num_tips = Rcpp::as<Rcpp::CharacterVector>(tree["tip.label"]).size();
  splittree::uvec node_names = splittree::Seq(static_cast<splittree::uint>(1), num_tips);
  
  vector<typename QuadraticPolynomialOU::LengthType> lengths(branches.n_rows);
  
  for(arma::uword i = 0; i < branches.n_rows; ++i) {
    lengths[i].length_ = t[i];
    lengths[i].regime_ = regimes[i] - 1;
  }
  
  if(threshold_SV <= 0) {
    ostringstream os;
    os<<"ERR:03812:PCMBaseCpp:Rcpp.cpp:CreateQuadraticPolynomialOU:: The argument threshold_SV should be positive real number.";
    throw invalid_argument(os.str());
  }
  
  
  if(threshold_Lambda_ij < 0) {
    ostringstream os;
    os<<"ERR:03813:PCMBaseCpp:Rcpp.cpp:CreateQuadraticPolynomialOU:: The argument threshold_Lambda_ij should be non-negative double.";
    throw invalid_argument(os.str());
  }
  
  typename QuadraticPolynomialOU::DataType data(
      node_names, X, Pc, internal_pc_full, 
      RModel, std::vector<std::string>(), 
      threshold_SV, threshold_skip_singular, skip_singular,
      threshold_Lambda_ij);
  
  return new QuadraticPolynomialOU(br_0, br_1, lengths, data);
}

RCPP_EXPOSED_CLASS_NODECL(QuadraticPolynomialOU::TraversalSpecificationType)
RCPP_EXPOSED_CLASS_NODECL(QuadraticPolynomialOU::AlgorithmType)
  
RCPP_MODULE(QuadraticPolynomialOU) {
  Rcpp::class_<QuadraticPolynomialOU::TreeType::Tree> ( "QuadraticPolynomialOU_Tree" )
  .property("num_nodes", &QuadraticPolynomialOU::TreeType::Tree::num_nodes )
  .property("num_tips", &QuadraticPolynomialOU::TreeType::Tree::num_tips )
  .method("FindNodeWithId", &QuadraticPolynomialOU::TreeType::Tree::FindNodeWithId )
  .method("FindIdOfNode", &QuadraticPolynomialOU::TreeType::Tree::FindIdOfNode )
  .method("FindIdOfParent", &QuadraticPolynomialOU::TreeType::Tree::FindIdOfParent )
  .method("OrderNodes", &QuadraticPolynomialOU::TreeType::Tree::OrderNodes )
  ;
  Rcpp::class_<QuadraticPolynomialOU::TreeType>( "QuadraticPolynomialOU_OrderedTree" )
    .derives<QuadraticPolynomialOU::TreeType::Tree> ( "QuadraticPolynomialOU_Tree" )
    .method("RangeIdPruneNode", &QuadraticPolynomialOU::TreeType::RangeIdPruneNode )
    .method("RangeIdVisitNode", &QuadraticPolynomialOU::TreeType::RangeIdVisitNode )
    .property("num_levels", &QuadraticPolynomialOU::TreeType::num_levels )
    .property("ranges_id_visit", &QuadraticPolynomialOU::TreeType::ranges_id_visit )
    .property("ranges_id_prune", &QuadraticPolynomialOU::TreeType::ranges_id_prune )
  ;
  Rcpp::class_<QuadraticPolynomialOU::AlgorithmType::ParentType>( "QuadraticPolynomialOU_TraversalAlgorithm" )
    .property( "VersionOPENMP", &QuadraticPolynomialOU::AlgorithmType::ParentType::VersionOPENMP )
    .property( "num_threads", &QuadraticPolynomialOU::AlgorithmType::num_threads )
  ;
  Rcpp::class_<QuadraticPolynomialOU::AlgorithmType> ( "QuadraticPolynomialOU_ParallelPruning" )
    .derives<QuadraticPolynomialOU::AlgorithmType::ParentType>( "QuadraticPolynomialOU_TraversalAlgorithm" )
    .method( "ModeAutoStep", &QuadraticPolynomialOU::AlgorithmType::ModeAutoStep )
    .property( "ModeAutoCurrent", &QuadraticPolynomialOU::AlgorithmType::ModeAutoCurrent )
    .property( "IsTuning", &QuadraticPolynomialOU::AlgorithmType::IsTuning )
    .property( "min_size_chunk_visit", &QuadraticPolynomialOU::AlgorithmType::min_size_chunk_visit )
    .property( "min_size_chunk_prune", &QuadraticPolynomialOU::AlgorithmType::min_size_chunk_prune )
    .property( "durations_tuning", &QuadraticPolynomialOU::AlgorithmType::durations_tuning )
    .property( "fastest_step_tuning", &QuadraticPolynomialOU::AlgorithmType::fastest_step_tuning )
  ;
  Rcpp::class_<QuadraticPolynomialOU::TraversalSpecificationType::BaseType> ( "QuadraticPolynomial_PrunignSpec" )
    .field( "A", &QuadraticPolynomialOU::TraversalSpecificationType::BaseType::A)
    .field( "b", &QuadraticPolynomialOU::TraversalSpecificationType::BaseType::b )
    .field( "C", &QuadraticPolynomialOU::TraversalSpecificationType::BaseType::C )
    .field( "d", &QuadraticPolynomialOU::TraversalSpecificationType::BaseType::d )
    .field( "E", &QuadraticPolynomialOU::TraversalSpecificationType::BaseType::E )
    .field( "f", &QuadraticPolynomialOU::TraversalSpecificationType::BaseType::f )
    .field( "L", &QuadraticPolynomialOU::TraversalSpecificationType::BaseType::L )
    .field( "m", &QuadraticPolynomialOU::TraversalSpecificationType::BaseType::m )
    .field( "r", &QuadraticPolynomialOU::TraversalSpecificationType::BaseType::r )
    .field( "V", &QuadraticPolynomialOU::TraversalSpecificationType::BaseType::V )
    .field( "V_1", &QuadraticPolynomialOU::TraversalSpecificationType::BaseType::V_1 )
  ;
  Rcpp::class_<QuadraticPolynomialOU::TraversalSpecificationType> ( "QuadraticPolynomialOU_PruningSpec" )
    .derives<QuadraticPolynomialOU::TraversalSpecificationType::BaseType>( "QuadraticPolynomial_PrunignSpec" )
  ;
  Rcpp::class_<QuadraticPolynomialOU>( "QuadraticPolynomialOU" )
    .factory<arma::mat const&, Rcpp::List const&, Rcpp::List const&>(&CreateQuadraticPolynomialOU)
    .method( "TraverseTree", &QuadraticPolynomialOU::TraverseTree )
    .property( "tree", &QuadraticPolynomialOU::tree )
    .property( "spec", &QuadraticPolynomialOU::spec )
    .property( "algorithm", &QuadraticPolynomialOU::algorithm )
  ;
}

QuadraticPolynomialJOU* CreateQuadraticPolynomialJOU(
    arma::mat const&X, 
    Rcpp::List const& tree, 
    Rcpp::List const& model,
    Rcpp::List const& metaInfo) {
  
  double threshold_SV = static_cast<double>(metaInfo["PCMBase.Threshold.SV"]);
  double threshold_skip_singular = static_cast<double>(metaInfo["PCMBase.Threshold.Skip.Singular"]);
  double threshold_Lambda_ij = static_cast<double>(metaInfo["PCMBase.Threshold.Lambda_ij"]);
  bool internal_pc_full = static_cast<int>(metaInfo["PCMBase.Internal.PC.Full"]);
  bool skip_singular = static_cast<int>(metaInfo["PCMBase.Skip.Singular"]);
  
  arma::imat Pc(X.n_rows, X.n_cols, arma::fill::ones);
  
  for(arma::uword i = 0; i < X.n_rows; ++i)
    for(arma::uword j = 0; j < X.n_cols; ++j) {
      Pc(i,j) = static_cast<arma::imat::value_type>(arma::is_finite(X(i,j)));
    }
    
    arma::umat branches = tree["edge"];
  splittree::uvec br_0 = arma::conv_to<splittree::uvec>::from(branches.col(0));
  splittree::uvec br_1 = arma::conv_to<splittree::uvec>::from(branches.col(1));
  splittree::vec t = Rcpp::as<splittree::vec>(tree["edge.length"]);
  
  using namespace std;
  
  uint RModel = Rcpp::as<uint>(metaInfo["RModel"]);
  vector<arma::uword> regimes = Rcpp::as<vector<arma::uword>>(metaInfo["r"]);
  vector<arma::u8> jumps = Rcpp::as<vector<arma::u8>>(metaInfo["xi"]);
  
  if(regimes.size() != branches.n_rows) {
    ostringstream os;
    os<<"ERR:03821:PCMBaseCpp:Rcpp.cpp:CreateQuadraticPolynomialJOU:: The slot r in metaInfo has different length ("<<regimes.size()<<
      ") than the number of edges ("<<branches.n_rows<<").";
    throw logic_error(os.str());
  }
  
  if(jumps.size() != branches.n_rows) {
    ostringstream os;
    os<<"ERR:03822:PCMBaseCpp:Rcpp.cpp:CreateQuadraticPolynomialJOU:: The slot jumps in trees has different length ("<<jumps.size()<<
      ") than the number of edges ("<<branches.n_rows<<").";
    throw logic_error(os.str());
  }
  
  splittree::uint num_tips = Rcpp::as<Rcpp::CharacterVector>(tree["tip.label"]).size();
  splittree::uvec node_names = splittree::Seq(static_cast<splittree::uint>(1), num_tips);
  
  vector<typename QuadraticPolynomialJOU::LengthType> lengths(branches.n_rows);
  
  for(arma::uword i = 0; i < branches.n_rows; ++i) {
    lengths[i].length_ = t[i];
    lengths[i].regime_ = regimes[i] - 1;
    lengths[i].jump_ = jumps[i];
  }
  
  if(threshold_SV <= 0) {
    ostringstream os;
    os<<"ERR:03823:PCMBaseCpp:Rcpp.cpp:CreateQuadraticPolynomialJOU:: The argument threshold_SV should be positive real number.";
    throw invalid_argument(os.str());
  }
  
  if(threshold_Lambda_ij < 0) {
    ostringstream os;
    os<<"ERR:03824:PCMBaseCpp:Rcpp.cpp:CreateQuadraticPolynomialJOU:: The argument threshold_Lambda_ij should be non-negative double.";
    throw invalid_argument(os.str());
  }
  
  typename QuadraticPolynomialJOU::DataType data(
      node_names, X, Pc, internal_pc_full, 
      RModel, std::vector<std::string>(), 
      threshold_SV, threshold_skip_singular, skip_singular,
      threshold_Lambda_ij);
  
  return new QuadraticPolynomialJOU(br_0, br_1, lengths, data);
}

RCPP_EXPOSED_CLASS_NODECL(QuadraticPolynomialJOU::TreeType)
RCPP_EXPOSED_CLASS_NODECL(QuadraticPolynomialJOU::TraversalSpecificationType)
RCPP_EXPOSED_CLASS_NODECL(QuadraticPolynomialJOU::AlgorithmType)
  
RCPP_MODULE(QuadraticPolynomialJOU) {
  Rcpp::class_<QuadraticPolynomialJOU::TreeType::Tree> ( "QuadraticPolynomialJOU_Tree" )
  .property("num_nodes", &QuadraticPolynomialJOU::TreeType::Tree::num_nodes )
  .property("num_tips", &QuadraticPolynomialJOU::TreeType::Tree::num_tips )
  .method("FindNodeWithId", &QuadraticPolynomialJOU::TreeType::Tree::FindNodeWithId )
  .method("FindIdOfNode", &QuadraticPolynomialJOU::TreeType::Tree::FindIdOfNode )
  .method("FindIdOfParent", &QuadraticPolynomialJOU::TreeType::Tree::FindIdOfParent )
  .method("OrderNodes", &QuadraticPolynomialJOU::TreeType::Tree::OrderNodes )
  ;
  Rcpp::class_<QuadraticPolynomialJOU::TreeType>( "QuadraticPolynomialJOU_OrderedTree" )
    .derives<QuadraticPolynomialJOU::TreeType::Tree> ( "QuadraticPolynomialJOU_Tree" )
    .method("RangeIdPruneNode", &QuadraticPolynomialJOU::TreeType::RangeIdPruneNode )
    .method("RangeIdVisitNode", &QuadraticPolynomialJOU::TreeType::RangeIdVisitNode )
    .property("num_levels", &QuadraticPolynomialJOU::TreeType::num_levels )
    .property("ranges_id_visit", &QuadraticPolynomialJOU::TreeType::ranges_id_visit )
    .property("ranges_id_prune", &QuadraticPolynomialJOU::TreeType::ranges_id_prune )
  ;
  Rcpp::class_<QuadraticPolynomialJOU::AlgorithmType::ParentType>( "QuadraticPolynomialJOU_TraversalAlgorithm" )
    .property( "VersionOPENMP", &QuadraticPolynomialJOU::AlgorithmType::ParentType::VersionOPENMP )
    .property( "num_threads", &QuadraticPolynomialJOU::AlgorithmType::num_threads )
  ;
  Rcpp::class_<QuadraticPolynomialJOU::AlgorithmType> ( "QuadraticPolynomialJOU_ParallelPruning" )
    .derives<QuadraticPolynomialJOU::AlgorithmType::ParentType>( "QuadraticPolynomialJOU_TraversalAlgorithm" )
    .method( "ModeAutoStep", &QuadraticPolynomialJOU::AlgorithmType::ModeAutoStep )
    .property( "ModeAutoCurrent", &QuadraticPolynomialJOU::AlgorithmType::ModeAutoCurrent )
    .property( "IsTuning", &QuadraticPolynomialJOU::AlgorithmType::IsTuning )
    .property( "min_size_chunk_visit", &QuadraticPolynomialJOU::AlgorithmType::min_size_chunk_visit )
    .property( "min_size_chunk_prune", &QuadraticPolynomialJOU::AlgorithmType::min_size_chunk_prune )
    .property( "durations_tuning", &QuadraticPolynomialJOU::AlgorithmType::durations_tuning )
    .property( "fastest_step_tuning", &QuadraticPolynomialJOU::AlgorithmType::fastest_step_tuning )
  ;
  Rcpp::class_<QuadraticPolynomialJOU::TraversalSpecificationType::BaseType> ( "QuadraticPolynomial_PrunignSpec" )
    .field( "A", &QuadraticPolynomialJOU::TraversalSpecificationType::BaseType::A)
    .field( "b", &QuadraticPolynomialJOU::TraversalSpecificationType::BaseType::b )
    .field( "C", &QuadraticPolynomialJOU::TraversalSpecificationType::BaseType::C )
    .field( "d", &QuadraticPolynomialJOU::TraversalSpecificationType::BaseType::d )
    .field( "E", &QuadraticPolynomialJOU::TraversalSpecificationType::BaseType::E )
    .field( "f", &QuadraticPolynomialJOU::TraversalSpecificationType::BaseType::f )
    .field( "L", &QuadraticPolynomialJOU::TraversalSpecificationType::BaseType::L )
    .field( "m", &QuadraticPolynomialJOU::TraversalSpecificationType::BaseType::m )
    .field( "r", &QuadraticPolynomialJOU::TraversalSpecificationType::BaseType::r )
    .field( "V", &QuadraticPolynomialJOU::TraversalSpecificationType::BaseType::V )
    .field( "V_1", &QuadraticPolynomialJOU::TraversalSpecificationType::BaseType::V_1 )
  ;
  Rcpp::class_<QuadraticPolynomialJOU::TraversalSpecificationType> ( "QuadraticPolynomialJOU_PruningSpec" )
    .derives<QuadraticPolynomialJOU::TraversalSpecificationType::BaseType>( "QuadraticPolynomial_PrunignSpec" )
  ;
  Rcpp::class_<QuadraticPolynomialJOU>( "QuadraticPolynomialJOU" )
    .factory<arma::mat const&, Rcpp::List const&, Rcpp::List const&>(&CreateQuadraticPolynomialJOU)
    .method( "TraverseTree", &QuadraticPolynomialJOU::TraverseTree )
    .property( "tree", &QuadraticPolynomialJOU::tree )
    .property( "spec", &QuadraticPolynomialJOU::spec )
    .property( "algorithm", &QuadraticPolynomialJOU::algorithm )
  ;
}

QuadraticPolynomialTwoSpeedOU* CreateQuadraticPolynomialTwoSpeedOU(
    arma::mat const& X,
    Rcpp::List const& tree,
    Rcpp::List const& model,
    Rcpp::List const& metaInfo) {
  
  double threshold_SV = static_cast<double>(metaInfo["PCMBase.Threshold.SV"]);
  double threshold_skip_singular = static_cast<double>(metaInfo["PCMBase.Threshold.Skip.Singular"]);
  double threshold_Lambda_ij = static_cast<double>(metaInfo["PCMBase.Threshold.Lambda_ij"]);
  bool internal_pc_full = static_cast<int>(metaInfo["PCMBase.Internal.PC.Full"]);
  bool skip_singular = static_cast<int>(metaInfo["PCMBase.Skip.Singular"]);
  
  arma::imat Pc(X.n_rows, X.n_cols, arma::fill::ones);
  
  for(arma::uword i = 0; i < X.n_rows; ++i)
    for(arma::uword j = 0; j < X.n_cols; ++j) {
      Pc(i,j) = static_cast<arma::imat::value_type>(arma::is_finite(X(i,j)));
    }
    
  arma::umat branches = tree["edge"];
  splittree::uvec br_0 = arma::conv_to<splittree::uvec>::from(branches.col(0));
  splittree::uvec br_1 = arma::conv_to<splittree::uvec>::from(branches.col(1));
  splittree::vec t = Rcpp::as<splittree::vec>(tree["edge.length"]);
  
  using namespace std;
  uint RModel = Rcpp::as<uint>(metaInfo["RModel"]);
  vector<arma::uword> regimes = Rcpp::as<vector<arma::uword>>(metaInfo["r"]);
  
  if(regimes.size() != branches.n_rows) {
    ostringstream os;
    os<<"ERR:03831:PCMBaseCpp:Rcpp.cpp:CreateQuadraticPolynomialTwoSpeedOU:: The slot r in metaInfo has different length ("<<regimes.size()<<
      ") than the number of edges ("<<branches.n_rows<<").";
    throw logic_error(os.str());
  }
  
  splittree::uint num_tips = Rcpp::as<Rcpp::CharacterVector>(tree["tip.label"]).size();
  splittree::uvec node_names = splittree::Seq(static_cast<splittree::uint>(1), num_tips);
  
  vector<typename QuadraticPolynomialTwoSpeedOU::LengthType> lengths(branches.n_rows);
  
  for(arma::uword i = 0; i < branches.n_rows; ++i) {
    lengths[i].length_ = t[i];
    lengths[i].regime_ = regimes[i] - 1;
  }
  
  if(threshold_SV <= 0) {
    ostringstream os;
    os<<"ERR:03832:PCMBaseCpp:Rcpp.cpp:CreateQuadraticPolynomialTwoSpeedOU:: The argument threshold_SV should be positive real number.";
    throw invalid_argument(os.str());
  }
  
  if(threshold_Lambda_ij < 0) {
    ostringstream os;
    os<<"ERR:03833:PCMBaseCpp:Rcpp.cpp:CreateQuadraticPolynomialTwoSpeedOU:: The argument threshold_Lambda_ij should be non-negative double.";
    throw invalid_argument(os.str());
  }
  
  typename QuadraticPolynomialTwoSpeedOU::DataType data(
      node_names, X, Pc, internal_pc_full, 
      RModel, std::vector<std::string>(), 
      threshold_SV, threshold_skip_singular, skip_singular,
      threshold_Lambda_ij);
  
  return new QuadraticPolynomialTwoSpeedOU(br_0, br_1, lengths, data);;
}

RCPP_EXPOSED_CLASS_NODECL(QuadraticPolynomialTwoSpeedOU::TraversalSpecificationType)
RCPP_EXPOSED_CLASS_NODECL(QuadraticPolynomialTwoSpeedOU::AlgorithmType)
  
RCPP_MODULE(QuadraticPolynomialTwoSpeedOU) {
  Rcpp::class_<QuadraticPolynomialTwoSpeedOU::TreeType::Tree> ( "QuadraticPolynomialTwoSpeedOU_Tree" )
  .property("num_nodes", &QuadraticPolynomialTwoSpeedOU::TreeType::Tree::num_nodes )
  .property("num_tips", &QuadraticPolynomialTwoSpeedOU::TreeType::Tree::num_tips )
  .method("FindNodeWithId", &QuadraticPolynomialTwoSpeedOU::TreeType::Tree::FindNodeWithId )
  .method("FindIdOfNode", &QuadraticPolynomialTwoSpeedOU::TreeType::Tree::FindIdOfNode )
  .method("FindIdOfParent", &QuadraticPolynomialTwoSpeedOU::TreeType::Tree::FindIdOfParent )
  .method("OrderNodes", &QuadraticPolynomialTwoSpeedOU::TreeType::Tree::OrderNodes )
  ;
  Rcpp::class_<QuadraticPolynomialTwoSpeedOU::TreeType>( "QuadraticPolynomialTwoSpeedOU_OrderedTree" )
    .derives<QuadraticPolynomialTwoSpeedOU::TreeType::Tree> ( "QuadraticPolynomialTwoSpeedOU_Tree" )
    .method("RangeIdPruneNode", &QuadraticPolynomialTwoSpeedOU::TreeType::RangeIdPruneNode )
    .method("RangeIdVisitNode", &QuadraticPolynomialTwoSpeedOU::TreeType::RangeIdVisitNode )
    .property("num_levels", &QuadraticPolynomialTwoSpeedOU::TreeType::num_levels )
    .property("ranges_id_visit", &QuadraticPolynomialTwoSpeedOU::TreeType::ranges_id_visit )
    .property("ranges_id_prune", &QuadraticPolynomialTwoSpeedOU::TreeType::ranges_id_prune )
  ;
  Rcpp::class_<QuadraticPolynomialTwoSpeedOU::AlgorithmType::ParentType>( "QuadraticPolynomialTwoSpeedOU_TraversalAlgorithm" )
    .property( "VersionOPENMP", &QuadraticPolynomialTwoSpeedOU::AlgorithmType::ParentType::VersionOPENMP )
    .property( "num_threads", &QuadraticPolynomialTwoSpeedOU::AlgorithmType::num_threads )
  ;
  Rcpp::class_<QuadraticPolynomialTwoSpeedOU::AlgorithmType> ( "QuadraticPolynomialTwoSpeedOU_ParallelPruning" )
    .derives<QuadraticPolynomialTwoSpeedOU::AlgorithmType::ParentType>( "QuadraticPolynomialTwoSpeedOU_TraversalAlgorithm" )
    .method( "ModeAutoStep", &QuadraticPolynomialTwoSpeedOU::AlgorithmType::ModeAutoStep )
    .property( "ModeAutoCurrent", &QuadraticPolynomialTwoSpeedOU::AlgorithmType::ModeAutoCurrent )
    .property( "IsTuning", &QuadraticPolynomialTwoSpeedOU::AlgorithmType::IsTuning )
    .property( "min_size_chunk_visit", &QuadraticPolynomialTwoSpeedOU::AlgorithmType::min_size_chunk_visit )
    .property( "min_size_chunk_prune", &QuadraticPolynomialTwoSpeedOU::AlgorithmType::min_size_chunk_prune )
    .property( "durations_tuning", &QuadraticPolynomialTwoSpeedOU::AlgorithmType::durations_tuning )
    .property( "fastest_step_tuning", &QuadraticPolynomialTwoSpeedOU::AlgorithmType::fastest_step_tuning )
  ;
  Rcpp::class_<QuadraticPolynomialTwoSpeedOU::TraversalSpecificationType::BaseType> ( "QuadraticPolynomial_PrunignSpec" )
    .field( "A", &QuadraticPolynomialTwoSpeedOU::TraversalSpecificationType::BaseType::A)
    .field( "b", &QuadraticPolynomialTwoSpeedOU::TraversalSpecificationType::BaseType::b )
    .field( "C", &QuadraticPolynomialTwoSpeedOU::TraversalSpecificationType::BaseType::C )
    .field( "d", &QuadraticPolynomialTwoSpeedOU::TraversalSpecificationType::BaseType::d )
    .field( "E", &QuadraticPolynomialTwoSpeedOU::TraversalSpecificationType::BaseType::E )
    .field( "f", &QuadraticPolynomialTwoSpeedOU::TraversalSpecificationType::BaseType::f )
    .field( "L", &QuadraticPolynomialTwoSpeedOU::TraversalSpecificationType::BaseType::L )
    .field( "m", &QuadraticPolynomialTwoSpeedOU::TraversalSpecificationType::BaseType::m )
    .field( "r", &QuadraticPolynomialTwoSpeedOU::TraversalSpecificationType::BaseType::r )
    .field( "V", &QuadraticPolynomialTwoSpeedOU::TraversalSpecificationType::BaseType::V )
    .field( "V_1", &QuadraticPolynomialTwoSpeedOU::TraversalSpecificationType::BaseType::V_1 )
  ;
  Rcpp::class_<QuadraticPolynomialTwoSpeedOU::TraversalSpecificationType> ( "QuadraticPolynomialTwoSpeedOU_PruningSpec" )
    .derives<QuadraticPolynomialTwoSpeedOU::TraversalSpecificationType::BaseType>( "QuadraticPolynomial_PrunignSpec" )
  ;
  Rcpp::class_<QuadraticPolynomialTwoSpeedOU>( "QuadraticPolynomialTwoSpeedOU" )
    .factory<arma::mat const&, Rcpp::List const&, Rcpp::List const&>(&CreateQuadraticPolynomialTwoSpeedOU)
    .method( "TraverseTree", &QuadraticPolynomialTwoSpeedOU::TraverseTree )
    .property( "tree", &QuadraticPolynomialTwoSpeedOU::tree )
    .property( "spec", &QuadraticPolynomialTwoSpeedOU::spec )
    .property( "algorithm", &QuadraticPolynomialTwoSpeedOU::algorithm )
  ;
}


QuadraticPolynomialMRG* CreateQuadraticPolynomialMRG(
    arma::mat const& X,
    Rcpp::List const& tree,
    Rcpp::List const& model,
    Rcpp::List const& metaInfo,
    std::vector<std::string> const& regimeModels) {
  
  double threshold_SV = static_cast<double>(metaInfo["PCMBase.Threshold.SV"]);
  double threshold_skip_singular = static_cast<double>(metaInfo["PCMBase.Threshold.Skip.Singular"]);
  double threshold_Lambda_ij = static_cast<double>(metaInfo["PCMBase.Threshold.Lambda_ij"]);
  bool internal_pc_full = static_cast<int>(metaInfo["PCMBase.Internal.PC.Full"]);
  bool skip_singular = static_cast<int>(metaInfo["PCMBase.Skip.Singular"]);
  
  // std::cout<<"skip_singular(1)"<<static_cast<bool>(metaInfo["PCMBase.Skip.Singular"])<<"\n";
  // std::cout<<"skip_singular(2)"<<static_cast<int>(metaInfo["PCMBase.Skip.Singular"])<<"\n";
  // std::cout<<"skip_singular(3)"<<skip_singular<<"\n";
  // 
  arma::imat Pc(X.n_rows, X.n_cols, arma::fill::ones);
  
  for(arma::uword i = 0; i < X.n_rows; ++i)
    for(arma::uword j = 0; j < X.n_cols; ++j) {
      Pc(i,j) = static_cast<arma::imat::value_type>(arma::is_finite(X(i,j)));
    }
    
  arma::umat branches = tree["edge"];
  splittree::uvec br_0 = arma::conv_to<splittree::uvec>::from(branches.col(0));
  splittree::uvec br_1 = arma::conv_to<splittree::uvec>::from(branches.col(1));
  splittree::vec t = Rcpp::as<splittree::vec>(tree["edge.length"]);
  
  using namespace std;
  uint RModel = Rcpp::as<uint>(metaInfo["RModel"]);
  vector<arma::uword> regimes = Rcpp::as<vector<arma::uword>>(metaInfo["r"]);
  vector<arma::u8> jumps = Rcpp::as<vector<arma::u8>>(metaInfo["xi"]);
  
  if(regimes.size() != branches.n_rows) {
    ostringstream os;
    os<<"ERR:03841:PCMBaseCpp:Rcpp.cpp:CreateQuadraticPolynomialMRG:: The slot r in metaInfo has different length ("<<regimes.size()<<
      ") than the number of edges ("<<branches.n_rows<<").";
    throw logic_error(os.str());
  }
  
  if(jumps.size() != branches.n_rows) {
    ostringstream os;
    os<<"ERR:03842:PCMBaseCpp:Rcpp.cpp:CreateQuadraticPolynomialJOU:: The slot jumps in trees has different length ("<<jumps.size()<<
      ") than the number of edges ("<<branches.n_rows<<").";
    throw logic_error(os.str());
  }
  
  splittree::uint num_tips = Rcpp::as<Rcpp::CharacterVector>(tree["tip.label"]).size();
  
  splittree::uvec node_names = splittree::Seq(static_cast<splittree::uint>(1), num_tips);
  
  vector<typename QuadraticPolynomialMRG::LengthType> lengths(branches.n_rows);
  
  for(arma::uword i = 0; i < branches.n_rows; ++i) {
    lengths[i].length_ = t[i];
    lengths[i].regime_ = regimes[i] - 1;
    lengths[i].jump_ = jumps[i];
  }
  
  
  if(threshold_SV <= 0) {
    ostringstream os;
    os<<"ERR:03843:PCMBaseCpp:Rcpp.cpp:CreateQuadraticPolynomialMRG:: The argument threshold_SV should be positive real number.";
    throw invalid_argument(os.str());
  }
  
  if(threshold_Lambda_ij < 0) {
    ostringstream os;
    os<<"ERR:03844:PCMBaseCpp:Rcpp.cpp:CreateQuadraticPolynomialMRG:: The argument threshold_Lambda_ij should be non-negative double.";
    throw invalid_argument(os.str());
  }
  
  typename QuadraticPolynomialMRG::DataType data(
      node_names, X, Pc, internal_pc_full, 
      RModel, 
      regimeModels,
      threshold_SV, threshold_skip_singular, skip_singular,
      threshold_Lambda_ij
      );
  
  return new QuadraticPolynomialMRG(br_0, br_1, lengths, data);
}

RCPP_EXPOSED_CLASS_NODECL(QuadraticPolynomialMRG::TraversalSpecificationType)
  RCPP_EXPOSED_CLASS_NODECL(QuadraticPolynomialMRG::AlgorithmType)
  
  RCPP_MODULE(QuadraticPolynomialMRG) {
    Rcpp::class_<QuadraticPolynomialMRG::TreeType::Tree> ( "QuadraticPolynomialMRG_Tree" )
    .property("num_nodes", &QuadraticPolynomialMRG::TreeType::Tree::num_nodes )
    .property("num_tips", &QuadraticPolynomialMRG::TreeType::Tree::num_tips )
    .method("FindNodeWithId", &QuadraticPolynomialMRG::TreeType::Tree::FindNodeWithId )
    .method("FindIdOfNode", &QuadraticPolynomialMRG::TreeType::Tree::FindIdOfNode )
    .method("FindIdOfParent", &QuadraticPolynomialMRG::TreeType::Tree::FindIdOfParent )
    .method("OrderNodes", &QuadraticPolynomialMRG::TreeType::Tree::OrderNodes )
    ;
    Rcpp::class_<QuadraticPolynomialMRG::TreeType>( "QuadraticPolynomialMRG_OrderedTree" )
      .derives<QuadraticPolynomialMRG::TreeType::Tree> ( "QuadraticPolynomialMRG_Tree" )
      .method("RangeIdPruneNode", &QuadraticPolynomialMRG::TreeType::RangeIdPruneNode )
      .method("RangeIdVisitNode", &QuadraticPolynomialMRG::TreeType::RangeIdVisitNode )
      .property("num_levels", &QuadraticPolynomialMRG::TreeType::num_levels )
      .property("ranges_id_visit", &QuadraticPolynomialMRG::TreeType::ranges_id_visit )
      .property("ranges_id_prune", &QuadraticPolynomialMRG::TreeType::ranges_id_prune )
    ;
    Rcpp::class_<QuadraticPolynomialMRG::AlgorithmType::ParentType>( "QuadraticPolynomialMRG_TraversalAlgorithm" )
      .property( "VersionOPENMP", &QuadraticPolynomialMRG::AlgorithmType::ParentType::VersionOPENMP )
      .property( "num_threads", &QuadraticPolynomialMRG::AlgorithmType::num_threads )
    ;
    Rcpp::class_<QuadraticPolynomialMRG::AlgorithmType> ( "QuadraticPolynomialMRG_ParallelPruning" )
      .derives<QuadraticPolynomialMRG::AlgorithmType::ParentType>( "QuadraticPolynomialMRG_TraversalAlgorithm" )
      .method( "ModeAutoStep", &QuadraticPolynomialMRG::AlgorithmType::ModeAutoStep )
      .property( "ModeAutoCurrent", &QuadraticPolynomialMRG::AlgorithmType::ModeAutoCurrent )
      .property( "IsTuning", &QuadraticPolynomialMRG::AlgorithmType::IsTuning )
      .property( "min_size_chunk_visit", &QuadraticPolynomialMRG::AlgorithmType::min_size_chunk_visit )
      .property( "min_size_chunk_prune", &QuadraticPolynomialMRG::AlgorithmType::min_size_chunk_prune )
      .property( "durations_tuning", &QuadraticPolynomialMRG::AlgorithmType::durations_tuning )
      .property( "fastest_step_tuning", &QuadraticPolynomialMRG::AlgorithmType::fastest_step_tuning )
    ;
    Rcpp::class_<QuadraticPolynomialMRG::TraversalSpecificationType::BaseType> ( "QuadraticPolynomial_PrunignSpec" )
      .field( "A", &QuadraticPolynomialMRG::TraversalSpecificationType::BaseType::A)
      .field( "b", &QuadraticPolynomialMRG::TraversalSpecificationType::BaseType::b )
      .field( "C", &QuadraticPolynomialMRG::TraversalSpecificationType::BaseType::C )
      .field( "d", &QuadraticPolynomialMRG::TraversalSpecificationType::BaseType::d )
      .field( "E", &QuadraticPolynomialMRG::TraversalSpecificationType::BaseType::E )
      .field( "f", &QuadraticPolynomialMRG::TraversalSpecificationType::BaseType::f )
      .field( "L", &QuadraticPolynomialMRG::TraversalSpecificationType::BaseType::L )
      .field( "m", &QuadraticPolynomialMRG::TraversalSpecificationType::BaseType::m )
      .field( "r", &QuadraticPolynomialMRG::TraversalSpecificationType::BaseType::r )
      .field( "V", &QuadraticPolynomialMRG::TraversalSpecificationType::BaseType::V )
      .field( "V_1", &QuadraticPolynomialMRG::TraversalSpecificationType::BaseType::V_1 )
    ;
    Rcpp::class_<QuadraticPolynomialMRG::TraversalSpecificationType> ( "QuadraticPolynomialMRG_PruningSpec" )
      .derives<QuadraticPolynomialMRG::TraversalSpecificationType::BaseType>( "QuadraticPolynomial_PrunignSpec" )
    ;
    Rcpp::class_<QuadraticPolynomialMRG>( "QuadraticPolynomialMRG" )
      .factory<arma::mat const&, Rcpp::List const&, Rcpp::List const&>(&CreateQuadraticPolynomialMRG)
      .method( "TraverseTree", &QuadraticPolynomialMRG::TraverseTree )
      .property( "tree", &QuadraticPolynomialMRG::tree )
      .property( "spec", &QuadraticPolynomialMRG::spec )
      .property( "algorithm", &QuadraticPolynomialMRG::algorithm )
    ;
  }
