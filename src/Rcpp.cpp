/*
 *  Rcpp.cpp
 *  PCMBaseCpp
 *
 * Copyright 2017,2018 Venelin Mitov
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

#include<vector>
#include<string>
#include<sstream>

#include "QuadraticPolynomialWhite.h"
#include "QuadraticPolynomialBM.h"
#include "QuadraticPolynomialOU.h"
#include "QuadraticPolynomialJOU.h"
#include "QuadraticPolynomialDOU.h"
#include "QuadraticPolynomialMixedGaussian.h"

// [[Rcpp::plugins("cpp11")]]
// [[Rcpp::plugins(openmp)]]
// [[Rcpp::depends(RcppArmadillo)]]

void R_init_PCMBaseCpp(DllInfo *info) {
   /* Register routines, allocate resources. */
   R_registerRoutines(info, NULL, NULL, NULL, NULL);
   R_useDynamicSymbols(info, TRUE);
}

void R_unload_PCMBaseCpp(DllInfo *info) {
   /* Release resources. */
}

using namespace PCMBaseCpp;
using namespace std;

SPLITT::Tree<uint, double>* CreatePCMBaseCppTree(Rcpp::List const& tree) {
  arma::umat branches = tree["edge"];
  SPLITT::uvec br_0 = arma::conv_to<SPLITT::uvec>::from(branches.col(0));
  SPLITT::uvec br_1 = arma::conv_to<SPLITT::uvec>::from(branches.col(1));
  SPLITT::vec t = Rcpp::as<SPLITT::vec>(tree["edge.length"]);
  return new SPLITT::Tree<uint, double>(br_0, br_1, t);
}

RCPP_MODULE(PCMBaseCpp__Tree) {
  Rcpp::class_<SPLITT::Tree<uint, double> > ( "PCMBaseCpp__Tree" )
  .factory<Rcpp::List const&>( &CreatePCMBaseCppTree )
  .property("num_nodes", &SPLITT::Tree<uint, double>::num_nodes )
  .property("num_tips", &SPLITT::Tree<uint, double>::num_tips )
  .method("LengthOfBranch", &SPLITT::Tree<uint, double>::LengthOfBranch )
  .method("FindNodeWithId", &SPLITT::Tree<uint, double>::FindNodeWithId )
  .method("FindIdOfNode", &SPLITT::Tree<uint, double>::FindIdOfNode )
  .method("FindIdOfParent", &SPLITT::Tree<uint, double>::FindIdOfParent )
  .method( "FindChildren", &SPLITT::Tree<uint, double>::FindChildren )
  .method("OrderNodes", &SPLITT::Tree<uint, double>::OrderNodes )
  ;
}

SPLITT::OrderedTree<uint, double>* CreatePCMBaseCppOrderedTree(Rcpp::List const& tree) {
  arma::umat branches = tree["edge"];
  SPLITT::uvec br_0 = arma::conv_to<SPLITT::uvec>::from(branches.col(0));
  SPLITT::uvec br_1 = arma::conv_to<SPLITT::uvec>::from(branches.col(1));
  SPLITT::vec t = Rcpp::as<SPLITT::vec>(tree["edge.length"]);
  return new SPLITT::OrderedTree<uint, double>(br_0, br_1, t);
}


RCPP_MODULE(PCMBaseCpp__OrderedTree) {
  Rcpp::class_<SPLITT::Tree<uint, double> > ( "PCMBaseCpp__Tree" )
  .factory<Rcpp::List const&>( &CreatePCMBaseCppTree )
  .property("num_nodes", &SPLITT::Tree<uint, double>::num_nodes )
  .property("num_tips", &SPLITT::Tree<uint, double>::num_tips )
  .method("LengthOfBranch", &SPLITT::Tree<uint, double>::LengthOfBranch )
  .method("FindNodeWithId", &SPLITT::Tree<uint, double>::FindNodeWithId )
  .method("FindIdOfNode", &SPLITT::Tree<uint, double>::FindIdOfNode )
  .method("FindIdOfParent", &SPLITT::Tree<uint, double>::FindIdOfParent )
  .method( "FindChildren", &SPLITT::Tree<uint, double>::FindChildren )
  .method("OrderNodes", &SPLITT::Tree<uint, double>::OrderNodes )
  ;
  Rcpp::class_<SPLITT::OrderedTree<uint, double> >( "PCMBaseCpp__OrderedTree" )
    .derives<SPLITT::Tree<uint, double> > ( "PCMBaseCpp__Tree" )
    .factory<Rcpp::List const&>( &CreatePCMBaseCppOrderedTree )
    .property("num_levels", &SPLITT::OrderedTree<uint, double>::num_levels )
    .property("num_parallel_ranges_prune", &SPLITT::OrderedTree<uint, double>::num_parallel_ranges_prune )
    .property("ranges_id_visit", &SPLITT::OrderedTree<uint, double>::ranges_id_visit )
    .property("ranges_id_prune", &SPLITT::OrderedTree<uint, double>::ranges_id_prune )
  ;
}

QuadraticPolynomialWhite* CreateQuadraticPolynomialWhite(
    arma::mat const&X, 
    Rcpp::List const& tree, 
    Rcpp::List const& model,
    Rcpp::List const& metaInfo) { 
  
  double threshold_SV = static_cast<double>(metaInfo["PCMBase.Threshold.SV"]);
  double threshold_EV = static_cast<double>(metaInfo["PCMBase.Threshold.EV"]);
  double threshold_skip_singular = static_cast<double>(metaInfo["PCMBase.Threshold.Skip.Singular"]);
  double threshold_Lambda_ij = static_cast<double>(metaInfo["PCMBase.Threshold.Lambda_ij"]);
  
  bool skip_singular = static_cast<int>(metaInfo["PCMBase.Skip.Singular"]);
  
  arma::mat SE(Rcpp::as<arma::mat>(metaInfo["SE"]));
    
  Rcpp::List pcListInt = Rcpp::as<Rcpp::List>(metaInfo["pcListInt"]);
  std::vector<arma::uvec> Pc(Rcpp::as<arma::uword>(metaInfo["M"]));
  for(arma::uword i = 0; i < Pc.size(); ++i) {
    Pc[i] = Rcpp::as<arma::uvec>(pcListInt[i]);
  }
    
  arma::umat branches = tree["edge"];
  SPLITT::uvec br_0 = arma::conv_to<SPLITT::uvec>::from(branches.col(0));
  SPLITT::uvec br_1 = arma::conv_to<SPLITT::uvec>::from(branches.col(1));
  SPLITT::vec t = Rcpp::as<SPLITT::vec>(tree["edge.length"]);
  
  using namespace std;
  
  uint RModel = Rcpp::as<uint>(metaInfo["RModel"]);
  
  vector<arma::uword> regimes = Rcpp::as<vector<arma::uword> >(metaInfo["r"]);
  
  if(regimes.size() != branches.n_rows) {
    ostringstream os;
    os<<"ERR:03801:PCMBaseCpp:Rcpp.cpp:CreateQuadraticPolynomialWhite:: The slot r in metaInfo has different length ("<<regimes.size()<<
      ") than the number of edges ("<<branches.n_rows<<").";
    throw logic_error(os.str());
  }
  
  SPLITT::uint num_tips = Rcpp::as<Rcpp::CharacterVector>(tree["tip.label"]).size();
  SPLITT::uvec tip_names = SPLITT::Seq(static_cast<SPLITT::uint>(1), num_tips);
  
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
  if(threshold_EV <= 0) {
    ostringstream os;
    os<<"ERR:03803:PCMBaseCpp:Rcpp.cpp:CreateQuadraticPolynomialMixedGaussian:: The argument threshold_EV should be positive real number.";
    throw invalid_argument(os.str());
  }
  typename QuadraticPolynomialWhite::DataType data(
      tip_names, X, SE, Pc, RModel, std::vector<std::string>(), 
      threshold_SV, threshold_EV, threshold_skip_singular, skip_singular,
      threshold_Lambda_ij);
  
  return new QuadraticPolynomialWhite(br_0, br_1, lengths, data);
}

  RCPP_EXPOSED_CLASS_NODECL(QuadraticPolynomialWhite::TreeType)
  RCPP_EXPOSED_CLASS_NODECL(QuadraticPolynomialWhite::TraversalSpecificationType)
  RCPP_EXPOSED_CLASS_NODECL(QuadraticPolynomialWhite::AlgorithmType)
  
  RCPP_MODULE(PCMBaseCpp__QuadraticPolynomialWhite) {
    Rcpp::class_<QuadraticPolynomialWhite::TreeType::Tree> ( "PCMBaseCpp__QuadraticPolynomialWhite_Tree" )
    .property("num_nodes", &QuadraticPolynomialWhite::TreeType::Tree::num_nodes )
    .property("num_tips", &QuadraticPolynomialWhite::TreeType::Tree::num_tips )
    .method("FindNodeWithId", &QuadraticPolynomialWhite::TreeType::Tree::FindNodeWithId )
    .method("FindIdOfNode", &QuadraticPolynomialWhite::TreeType::Tree::FindIdOfNode )
    .method("FindIdOfParent", &QuadraticPolynomialWhite::TreeType::Tree::FindIdOfParent )
    .method("OrderNodes", &QuadraticPolynomialWhite::TreeType::Tree::OrderNodes )
    ;
    Rcpp::class_<QuadraticPolynomialWhite::TreeType>( "PCMBaseCpp__QuadraticPolynomialWhite_OrderedTree" )
      .derives<QuadraticPolynomialWhite::TreeType::Tree> ( "PCMBaseCpp__QuadraticPolynomialWhite_Tree" )
      .method("RangeIdPruneNode", &QuadraticPolynomialWhite::TreeType::RangeIdPruneNode )
      .method("RangeIdVisitNode", &QuadraticPolynomialWhite::TreeType::RangeIdVisitNode )
      .property("num_levels", &QuadraticPolynomialWhite::TreeType::num_levels )
      .property("ranges_id_visit", &QuadraticPolynomialWhite::TreeType::ranges_id_visit )
      .property("ranges_id_prune", &QuadraticPolynomialWhite::TreeType::ranges_id_prune )
    ;
    Rcpp::class_<QuadraticPolynomialWhite::AlgorithmType::ParentType>( "PCMBaseCpp__QuadraticPolynomialWhite_TraversalAlgorithm" )
      .property( "VersionOPENMP", &QuadraticPolynomialWhite::AlgorithmType::ParentType::VersionOPENMP )
      .property( "NumOmpThreads", &QuadraticPolynomialWhite::AlgorithmType::NumOmpThreads )
    ;
    Rcpp::class_<QuadraticPolynomialWhite::AlgorithmType> ( "PCMBaseCpp__QuadraticPolynomialWhite_ParallelPruning" )
      .derives<QuadraticPolynomialWhite::AlgorithmType::ParentType>( "PCMBaseCpp__QuadraticPolynomialWhite_TraversalAlgorithm" )
      .method( "ModeAutoStep", &QuadraticPolynomialWhite::AlgorithmType::ModeAutoStep )
      .property( "ModeAutoCurrent", &QuadraticPolynomialWhite::AlgorithmType::ModeAutoCurrent )
      .property( "IsTuning", &QuadraticPolynomialWhite::AlgorithmType::IsTuning )
      .property( "min_size_chunk_visit", &QuadraticPolynomialWhite::AlgorithmType::min_size_chunk_visit )
      .property( "min_size_chunk_prune", &QuadraticPolynomialWhite::AlgorithmType::min_size_chunk_prune )
      .property( "durations_tuning", &QuadraticPolynomialWhite::AlgorithmType::durations_tuning )
      .property( "fastest_step_tuning", &QuadraticPolynomialWhite::AlgorithmType::fastest_step_tuning )
    ;
    Rcpp::class_<QuadraticPolynomialWhite::TraversalSpecificationType::BaseType> ( "PCMBaseCpp__QuadraticPolynomial_PrunignSpec" )
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
    Rcpp::class_<QuadraticPolynomialWhite::TraversalSpecificationType> ( "PCMBaseCpp__QuadraticPolynomialWhite_PruningSpec" )
      .derives<QuadraticPolynomialWhite::TraversalSpecificationType::BaseType>( "PCMBaseCpp__QuadraticPolynomial_PrunignSpec" )
    ;
    Rcpp::class_<QuadraticPolynomialWhite>( "PCMBaseCpp__QuadraticPolynomialWhite" )
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
  double threshold_EV = static_cast<double>(metaInfo["PCMBase.Threshold.EV"]);
  double threshold_skip_singular = static_cast<double>(metaInfo["PCMBase.Threshold.Skip.Singular"]);
  double threshold_Lambda_ij = static_cast<double>(metaInfo["PCMBase.Threshold.Lambda_ij"]);
  
  bool skip_singular = static_cast<int>(metaInfo["PCMBase.Skip.Singular"]);
  
  arma::mat SE(Rcpp::as<arma::mat>(metaInfo["SE"]));
  
  Rcpp::List pcListInt = Rcpp::as<Rcpp::List>(metaInfo["pcListInt"]);
  std::vector<arma::uvec> Pc(Rcpp::as<arma::uword>(metaInfo["M"]));
  for(arma::uword i = 0; i < Pc.size(); ++i) {
    Pc[i] = Rcpp::as<arma::uvec>(pcListInt[i]);
  }
    
  arma::umat branches = tree["edge"];
  SPLITT::uvec br_0 = arma::conv_to<SPLITT::uvec>::from(branches.col(0));
  SPLITT::uvec br_1 = arma::conv_to<SPLITT::uvec>::from(branches.col(1));
  SPLITT::vec t = Rcpp::as<SPLITT::vec>(tree["edge.length"]);
  
  using namespace std;
  
  uint RModel = Rcpp::as<uint>(metaInfo["RModel"]);
  
  vector<arma::uword> regimes = Rcpp::as<vector<arma::uword> >(metaInfo["r"]);
  
  if(regimes.size() != branches.n_rows) {
    ostringstream os;
    os<<"ERR:03801:PCMBaseCpp:Rcpp.cpp:CreateQuadraticPolynomialBM:: The slot r in metaInfo has different length ("<<regimes.size()<<
      ") than the number of edges ("<<branches.n_rows<<").";
    throw logic_error(os.str());
  }
  
  SPLITT::uint num_tips = Rcpp::as<Rcpp::CharacterVector>(tree["tip.label"]).size();
  SPLITT::uvec tip_names = SPLITT::Seq(static_cast<SPLITT::uint>(1), num_tips);
  
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
  if(threshold_EV <= 0) {
    ostringstream os;
    os<<"ERR:03803:PCMBaseCpp:Rcpp.cpp:CreateQuadraticPolynomialMixedGaussian:: The argument threshold_EV should be positive real number.";
    throw invalid_argument(os.str());
  }
  typename QuadraticPolynomialBM::DataType data(
      tip_names, X, SE, Pc, RModel, std::vector<std::string>(), 
      threshold_SV, threshold_EV, threshold_skip_singular, skip_singular,
      threshold_Lambda_ij);
  
  return new QuadraticPolynomialBM(br_0, br_1, lengths, data);
}

//RCPP_EXPOSED_CLASS_NODECL(QuadraticPolynomialBM::TreeType)
RCPP_EXPOSED_CLASS_NODECL(QuadraticPolynomialBM::TraversalSpecificationType)
RCPP_EXPOSED_CLASS_NODECL(QuadraticPolynomialBM::AlgorithmType)
  
RCPP_MODULE(PCMBaseCpp__QuadraticPolynomialBM) {
  Rcpp::class_<QuadraticPolynomialBM::TreeType::Tree> ( "PCMBaseCpp__QuadraticPolynomialBM_Tree" )
  .property("num_nodes", &QuadraticPolynomialBM::TreeType::Tree::num_nodes )
  .property("num_tips", &QuadraticPolynomialBM::TreeType::Tree::num_tips )
  .method("FindNodeWithId", &QuadraticPolynomialBM::TreeType::Tree::FindNodeWithId )
  .method("FindIdOfNode", &QuadraticPolynomialBM::TreeType::Tree::FindIdOfNode )
  .method("FindIdOfParent", &QuadraticPolynomialBM::TreeType::Tree::FindIdOfParent )
  .method("OrderNodes", &QuadraticPolynomialBM::TreeType::Tree::OrderNodes )
  ;
  Rcpp::class_<QuadraticPolynomialBM::TreeType>( "PCMBaseCpp__QuadraticPolynomialBM_OrderedTree" )
    .derives<QuadraticPolynomialBM::TreeType::Tree> ( "PCMBaseCpp__QuadraticPolynomialBM_Tree" )
    .method("RangeIdPruneNode", &QuadraticPolynomialBM::TreeType::RangeIdPruneNode )
    .method("RangeIdVisitNode", &QuadraticPolynomialBM::TreeType::RangeIdVisitNode )
    .property("num_levels", &QuadraticPolynomialBM::TreeType::num_levels )
    .property("ranges_id_visit", &QuadraticPolynomialBM::TreeType::ranges_id_visit )
    .property("ranges_id_prune", &QuadraticPolynomialBM::TreeType::ranges_id_prune )
  ;
  Rcpp::class_<QuadraticPolynomialBM::AlgorithmType::ParentType>( "PCMBaseCpp__QuadraticPolynomialBM_TraversalAlgorithm" )
    .property( "VersionOPENMP", &QuadraticPolynomialBM::AlgorithmType::ParentType::VersionOPENMP )
    .property( "NumOmpThreads", &QuadraticPolynomialBM::AlgorithmType::NumOmpThreads )
  ;
  Rcpp::class_<QuadraticPolynomialBM::AlgorithmType> ( "PCMBaseCpp__QuadraticPolynomialBM_ParallelPruning" )
    .derives<QuadraticPolynomialBM::AlgorithmType::ParentType>( "PCMBaseCpp__QuadraticPolynomialBM_TraversalAlgorithm" )
    .method( "ModeAutoStep", &QuadraticPolynomialBM::AlgorithmType::ModeAutoStep )
    .property( "ModeAutoCurrent", &QuadraticPolynomialBM::AlgorithmType::ModeAutoCurrent )
    .property( "IsTuning", &QuadraticPolynomialBM::AlgorithmType::IsTuning )
    .property( "min_size_chunk_visit", &QuadraticPolynomialBM::AlgorithmType::min_size_chunk_visit )
    .property( "min_size_chunk_prune", &QuadraticPolynomialBM::AlgorithmType::min_size_chunk_prune )
    .property( "durations_tuning", &QuadraticPolynomialBM::AlgorithmType::durations_tuning )
    .property( "fastest_step_tuning", &QuadraticPolynomialBM::AlgorithmType::fastest_step_tuning )
  ;
  Rcpp::class_<QuadraticPolynomialBM::TraversalSpecificationType::BaseType> ( "PCMBaseCpp__QuadraticPolynomial_PrunignSpec" )
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
  Rcpp::class_<QuadraticPolynomialBM::TraversalSpecificationType> ( "PCMBaseCpp__QuadraticPolynomialBM_PruningSpec" )
    .derives<QuadraticPolynomialBM::TraversalSpecificationType::BaseType>( "PCMBaseCpp__QuadraticPolynomial_PrunignSpec" )
  ;
  Rcpp::class_<QuadraticPolynomialBM>( "PCMBaseCpp__QuadraticPolynomialBM" )
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
  double threshold_EV = static_cast<double>(metaInfo["PCMBase.Threshold.EV"]);
  double threshold_skip_singular = static_cast<double>(metaInfo["PCMBase.Threshold.Skip.Singular"]);
  double threshold_Lambda_ij = static_cast<double>(metaInfo["PCMBase.Threshold.Lambda_ij"]);
  
  bool skip_singular = static_cast<int>(metaInfo["PCMBase.Skip.Singular"]);
  
  arma::mat SE(Rcpp::as<arma::mat>(metaInfo["SE"]));
  
  Rcpp::List pcListInt = Rcpp::as<Rcpp::List>(metaInfo["pcListInt"]);
  std::vector<arma::uvec> Pc(Rcpp::as<arma::uword>(metaInfo["M"]));
  for(arma::uword i = 0; i < Pc.size(); ++i) {
    Pc[i] = Rcpp::as<arma::uvec>(pcListInt[i]);
  }
  
  arma::umat branches = tree["edge"];
  SPLITT::uvec br_0 = arma::conv_to<SPLITT::uvec>::from(branches.col(0));
  SPLITT::uvec br_1 = arma::conv_to<SPLITT::uvec>::from(branches.col(1));
  SPLITT::vec t = Rcpp::as<SPLITT::vec>(tree["edge.length"]);
  
  using namespace std;

  uint RModel = Rcpp::as<uint>(metaInfo["RModel"]);
  
  vector<arma::uword> regimes = Rcpp::as<vector<arma::uword> >(metaInfo["r"]);
  

  if(regimes.size() != branches.n_rows) {
    ostringstream os;
    os<<"ERR:03811:PCMBaseCpp:Rcpp.cpp:CreateQuadraticPolynomialOU:: The slot r in metaInfo has different length ("<<regimes.size()<<
      ") than the number of edges ("<<branches.n_rows<<").";
    throw logic_error(os.str());
  }

  SPLITT::uint num_tips = Rcpp::as<Rcpp::CharacterVector>(tree["tip.label"]).size();
  SPLITT::uvec tip_names = SPLITT::Seq(static_cast<SPLITT::uint>(1), num_tips);
  
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
  if(threshold_EV <= 0) {
    ostringstream os;
    os<<"ERR:03814:PCMBaseCpp:Rcpp.cpp:CreateQuadraticPolynomialMixedGaussian:: The argument threshold_EV should be positive real number.";
    throw invalid_argument(os.str());
  }
  
  if(threshold_Lambda_ij < 0) {
    ostringstream os;
    os<<"ERR:03813:PCMBaseCpp:Rcpp.cpp:CreateQuadraticPolynomialOU:: The argument threshold_Lambda_ij should be non-negative double.";
    throw invalid_argument(os.str());
  }
  
  typename QuadraticPolynomialOU::DataType data(
      tip_names, X, SE, Pc, 
      RModel, std::vector<std::string>(), 
      threshold_SV, threshold_EV, threshold_skip_singular, skip_singular,
      threshold_Lambda_ij);
  
  return new QuadraticPolynomialOU(br_0, br_1, lengths, data);
}

RCPP_EXPOSED_CLASS_NODECL(QuadraticPolynomialOU::TraversalSpecificationType)
RCPP_EXPOSED_CLASS_NODECL(QuadraticPolynomialOU::AlgorithmType)
  
RCPP_MODULE(PCMBaseCpp__QuadraticPolynomialOU) {
  Rcpp::class_<QuadraticPolynomialOU::TreeType::Tree> ( "PCMBaseCpp__QuadraticPolynomialOU_Tree" )
  .property("num_nodes", &QuadraticPolynomialOU::TreeType::Tree::num_nodes )
  .property("num_tips", &QuadraticPolynomialOU::TreeType::Tree::num_tips )
  .method("FindNodeWithId", &QuadraticPolynomialOU::TreeType::Tree::FindNodeWithId )
  .method("FindIdOfNode", &QuadraticPolynomialOU::TreeType::Tree::FindIdOfNode )
  .method("FindIdOfParent", &QuadraticPolynomialOU::TreeType::Tree::FindIdOfParent )
  .method("OrderNodes", &QuadraticPolynomialOU::TreeType::Tree::OrderNodes )
  ;
  Rcpp::class_<QuadraticPolynomialOU::TreeType>( "PCMBaseCpp__QuadraticPolynomialOU_OrderedTree" )
    .derives<QuadraticPolynomialOU::TreeType::Tree> ( "PCMBaseCpp__QuadraticPolynomialOU_Tree" )
    .method("RangeIdPruneNode", &QuadraticPolynomialOU::TreeType::RangeIdPruneNode )
    .method("RangeIdVisitNode", &QuadraticPolynomialOU::TreeType::RangeIdVisitNode )
    .property("num_levels", &QuadraticPolynomialOU::TreeType::num_levels )
    .property("ranges_id_visit", &QuadraticPolynomialOU::TreeType::ranges_id_visit )
    .property("ranges_id_prune", &QuadraticPolynomialOU::TreeType::ranges_id_prune )
  ;
  Rcpp::class_<QuadraticPolynomialOU::AlgorithmType::ParentType>( "PCMBaseCpp__QuadraticPolynomialOU_TraversalAlgorithm" )
    .property( "VersionOPENMP", &QuadraticPolynomialOU::AlgorithmType::ParentType::VersionOPENMP )
    .property( "NumOmpThreads", &QuadraticPolynomialOU::AlgorithmType::NumOmpThreads )
  ;
  Rcpp::class_<QuadraticPolynomialOU::AlgorithmType> ( "PCMBaseCpp__QuadraticPolynomialOU_ParallelPruning" )
    .derives<QuadraticPolynomialOU::AlgorithmType::ParentType>( "PCMBaseCpp__QuadraticPolynomialOU_TraversalAlgorithm" )
    .method( "ModeAutoStep", &QuadraticPolynomialOU::AlgorithmType::ModeAutoStep )
    .property( "ModeAutoCurrent", &QuadraticPolynomialOU::AlgorithmType::ModeAutoCurrent )
    .property( "IsTuning", &QuadraticPolynomialOU::AlgorithmType::IsTuning )
    .property( "min_size_chunk_visit", &QuadraticPolynomialOU::AlgorithmType::min_size_chunk_visit )
    .property( "min_size_chunk_prune", &QuadraticPolynomialOU::AlgorithmType::min_size_chunk_prune )
    .property( "durations_tuning", &QuadraticPolynomialOU::AlgorithmType::durations_tuning )
    .property( "fastest_step_tuning", &QuadraticPolynomialOU::AlgorithmType::fastest_step_tuning )
  ;
  Rcpp::class_<QuadraticPolynomialOU::TraversalSpecificationType::BaseType> ( "PCMBaseCpp__QuadraticPolynomial_PrunignSpec" )
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
  Rcpp::class_<QuadraticPolynomialOU::TraversalSpecificationType> ( "PCMBaseCpp__QuadraticPolynomialOU_PruningSpec" )
    .derives<QuadraticPolynomialOU::TraversalSpecificationType::BaseType>( "PCMBaseCpp__QuadraticPolynomial_PrunignSpec" )
  ;
  Rcpp::class_<QuadraticPolynomialOU>( "PCMBaseCpp__QuadraticPolynomialOU" )
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
  double threshold_EV = static_cast<double>(metaInfo["PCMBase.Threshold.EV"]);
  double threshold_skip_singular = static_cast<double>(metaInfo["PCMBase.Threshold.Skip.Singular"]);
  double threshold_Lambda_ij = static_cast<double>(metaInfo["PCMBase.Threshold.Lambda_ij"]);
  
  bool skip_singular = static_cast<int>(metaInfo["PCMBase.Skip.Singular"]);
 
  arma::mat SE(Rcpp::as<arma::mat>(metaInfo["SE"]));
  Rcpp::List pcListInt = Rcpp::as<Rcpp::List>(metaInfo["pcListInt"]);
  std::vector<arma::uvec> Pc(Rcpp::as<arma::uword>(metaInfo["M"]));
  for(arma::uword i = 0; i < Pc.size(); ++i) {
    Pc[i] = Rcpp::as<arma::uvec>(pcListInt[i]);
  }
 
  arma::umat branches = tree["edge"];
  SPLITT::uvec br_0 = arma::conv_to<SPLITT::uvec>::from(branches.col(0));
  SPLITT::uvec br_1 = arma::conv_to<SPLITT::uvec>::from(branches.col(1));
  SPLITT::vec t = Rcpp::as<SPLITT::vec>(tree["edge.length"]);
  
  using namespace std;
  
  uint RModel = Rcpp::as<uint>(metaInfo["RModel"]);
  vector<arma::uword> regimes = Rcpp::as<vector<arma::uword> >(metaInfo["r"]);
  vector<arma::u8> jumps = Rcpp::as<vector<arma::u8> >(metaInfo["xi"]);
  
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
  
  SPLITT::uint num_tips = Rcpp::as<Rcpp::CharacterVector>(tree["tip.label"]).size();
  SPLITT::uvec tip_names = SPLITT::Seq(static_cast<SPLITT::uint>(1), num_tips);
  
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
  if(threshold_EV <= 0) {
    ostringstream os;
    os<<"ERR:03824:PCMBaseCpp:Rcpp.cpp:CreateQuadraticPolynomialMixedGaussian:: The argument threshold_EV should be positive real number.";
    throw invalid_argument(os.str());
  }
  if(threshold_Lambda_ij < 0) {
    ostringstream os;
    os<<"ERR:03825:PCMBaseCpp:Rcpp.cpp:CreateQuadraticPolynomialJOU:: The argument threshold_Lambda_ij should be non-negative double.";
    throw invalid_argument(os.str());
  }
  
  typename QuadraticPolynomialJOU::DataType data(
      tip_names, X, SE, Pc, 
      RModel, std::vector<std::string>(), 
      threshold_SV, threshold_EV, threshold_skip_singular, skip_singular,
      threshold_Lambda_ij);
  
  return new QuadraticPolynomialJOU(br_0, br_1, lengths, data);
}

RCPP_EXPOSED_CLASS_NODECL(QuadraticPolynomialJOU::TreeType)
RCPP_EXPOSED_CLASS_NODECL(QuadraticPolynomialJOU::TraversalSpecificationType)
RCPP_EXPOSED_CLASS_NODECL(QuadraticPolynomialJOU::AlgorithmType)
  
RCPP_MODULE(PCMBaseCpp__QuadraticPolynomialJOU) {
  Rcpp::class_<QuadraticPolynomialJOU::TreeType::Tree> ( "PCMBaseCpp__QuadraticPolynomialJOU_Tree" )
  .property("num_nodes", &QuadraticPolynomialJOU::TreeType::Tree::num_nodes )
  .property("num_tips", &QuadraticPolynomialJOU::TreeType::Tree::num_tips )
  .method("FindNodeWithId", &QuadraticPolynomialJOU::TreeType::Tree::FindNodeWithId )
  .method("FindIdOfNode", &QuadraticPolynomialJOU::TreeType::Tree::FindIdOfNode )
  .method("FindIdOfParent", &QuadraticPolynomialJOU::TreeType::Tree::FindIdOfParent )
  .method("OrderNodes", &QuadraticPolynomialJOU::TreeType::Tree::OrderNodes )
  ;
  Rcpp::class_<QuadraticPolynomialJOU::TreeType>( "PCMBaseCpp__QuadraticPolynomialJOU_OrderedTree" )
    .derives<QuadraticPolynomialJOU::TreeType::Tree> ( "PCMBaseCpp__QuadraticPolynomialJOU_Tree" )
    .method("RangeIdPruneNode", &QuadraticPolynomialJOU::TreeType::RangeIdPruneNode )
    .method("RangeIdVisitNode", &QuadraticPolynomialJOU::TreeType::RangeIdVisitNode )
    .property("num_levels", &QuadraticPolynomialJOU::TreeType::num_levels )
    .property("ranges_id_visit", &QuadraticPolynomialJOU::TreeType::ranges_id_visit )
    .property("ranges_id_prune", &QuadraticPolynomialJOU::TreeType::ranges_id_prune )
  ;
  Rcpp::class_<QuadraticPolynomialJOU::AlgorithmType::ParentType>( "PCMBaseCpp__QuadraticPolynomialJOU_TraversalAlgorithm" )
    .property( "VersionOPENMP", &QuadraticPolynomialJOU::AlgorithmType::ParentType::VersionOPENMP )
    .property( "NumOmpThreads", &QuadraticPolynomialJOU::AlgorithmType::NumOmpThreads )
  ;
  Rcpp::class_<QuadraticPolynomialJOU::AlgorithmType> ( "PCMBaseCpp__QuadraticPolynomialJOU_ParallelPruning" )
    .derives<QuadraticPolynomialJOU::AlgorithmType::ParentType>( "PCMBaseCpp__QuadraticPolynomialJOU_TraversalAlgorithm" )
    .method( "ModeAutoStep", &QuadraticPolynomialJOU::AlgorithmType::ModeAutoStep )
    .property( "ModeAutoCurrent", &QuadraticPolynomialJOU::AlgorithmType::ModeAutoCurrent )
    .property( "IsTuning", &QuadraticPolynomialJOU::AlgorithmType::IsTuning )
    .property( "min_size_chunk_visit", &QuadraticPolynomialJOU::AlgorithmType::min_size_chunk_visit )
    .property( "min_size_chunk_prune", &QuadraticPolynomialJOU::AlgorithmType::min_size_chunk_prune )
    .property( "durations_tuning", &QuadraticPolynomialJOU::AlgorithmType::durations_tuning )
    .property( "fastest_step_tuning", &QuadraticPolynomialJOU::AlgorithmType::fastest_step_tuning )
  ;
  Rcpp::class_<QuadraticPolynomialJOU::TraversalSpecificationType::BaseType> ( "PCMBaseCpp__QuadraticPolynomial_PrunignSpec" )
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
  Rcpp::class_<QuadraticPolynomialJOU::TraversalSpecificationType> ( "PCMBaseCpp__QuadraticPolynomialJOU_PruningSpec" )
    .derives<QuadraticPolynomialJOU::TraversalSpecificationType::BaseType>( "PCMBaseCpp__QuadraticPolynomial_PrunignSpec" )
  ;
  Rcpp::class_<QuadraticPolynomialJOU>( "PCMBaseCpp__QuadraticPolynomialJOU" )
    .factory<arma::mat const&, Rcpp::List const&, Rcpp::List const&>(&CreateQuadraticPolynomialJOU)
    .method( "TraverseTree", &QuadraticPolynomialJOU::TraverseTree )
    .property( "tree", &QuadraticPolynomialJOU::tree )
    .property( "spec", &QuadraticPolynomialJOU::spec )
    .property( "algorithm", &QuadraticPolynomialJOU::algorithm )
  ;
}

QuadraticPolynomialDOU* CreateQuadraticPolynomialDOU(
    arma::mat const& X,
    Rcpp::List const& tree,
    Rcpp::List const& model,
    Rcpp::List const& metaInfo) {
  
  double threshold_SV = static_cast<double>(metaInfo["PCMBase.Threshold.SV"]);
  double threshold_EV = static_cast<double>(metaInfo["PCMBase.Threshold.EV"]);
  double threshold_skip_singular = static_cast<double>(metaInfo["PCMBase.Threshold.Skip.Singular"]);
  double threshold_Lambda_ij = static_cast<double>(metaInfo["PCMBase.Threshold.Lambda_ij"]);
  
  bool skip_singular = static_cast<int>(metaInfo["PCMBase.Skip.Singular"]);

  arma::mat SE(Rcpp::as<arma::mat>(metaInfo["SE"]));
  
  Rcpp::List pcListInt = Rcpp::as<Rcpp::List>(metaInfo["pcListInt"]);
  std::vector<arma::uvec> Pc(Rcpp::as<arma::uword>(metaInfo["M"]));
  for(arma::uword i = 0; i < Pc.size(); ++i) {
    Pc[i] = Rcpp::as<arma::uvec>(pcListInt[i]);
  }
  
  arma::umat branches = tree["edge"];
  SPLITT::uvec br_0 = arma::conv_to<SPLITT::uvec>::from(branches.col(0));
  SPLITT::uvec br_1 = arma::conv_to<SPLITT::uvec>::from(branches.col(1));
  SPLITT::vec t = Rcpp::as<SPLITT::vec>(tree["edge.length"]);
  
  using namespace std;
  uint RModel = Rcpp::as<uint>(metaInfo["RModel"]);
  vector<arma::uword> regimes = Rcpp::as<vector<arma::uword> >(metaInfo["r"]);
  
  if(regimes.size() != branches.n_rows) {
    ostringstream os;
    os<<"ERR:03831:PCMBaseCpp:Rcpp.cpp:CreateQuadraticPolynomialDOU:: The slot r in metaInfo has different length ("<<regimes.size()<<
      ") than the number of edges ("<<branches.n_rows<<").";
    throw logic_error(os.str());
  }
  
  SPLITT::uint num_tips = Rcpp::as<Rcpp::CharacterVector>(tree["tip.label"]).size();
  SPLITT::uvec tip_names = SPLITT::Seq(static_cast<SPLITT::uint>(1), num_tips);
  
  vector<typename QuadraticPolynomialDOU::LengthType> lengths(branches.n_rows);
  
  for(arma::uword i = 0; i < branches.n_rows; ++i) {
    lengths[i].length_ = t[i];
    lengths[i].regime_ = regimes[i] - 1;
  }
  
  if(threshold_SV <= 0) {
    ostringstream os;
    os<<"ERR:03832:PCMBaseCpp:Rcpp.cpp:CreateQuadraticPolynomialDOU:: The argument threshold_SV should be positive real number.";
    throw invalid_argument(os.str());
  }
  if(threshold_EV <= 0) {
    ostringstream os;
    os<<"ERR:03834:PCMBaseCpp:Rcpp.cpp:CreateQuadraticPolynomialMixedGaussian:: The argument threshold_EV should be positive real number.";
    throw invalid_argument(os.str());
  }
  if(threshold_Lambda_ij < 0) {
    ostringstream os;
    os<<"ERR:03833:PCMBaseCpp:Rcpp.cpp:CreateQuadraticPolynomialDOU:: The argument threshold_Lambda_ij should be non-negative double.";
    throw invalid_argument(os.str());
  }
  
  typename QuadraticPolynomialDOU::DataType data(
      tip_names, X, SE, Pc, 
      RModel, std::vector<std::string>(), 
      threshold_SV, threshold_EV, threshold_skip_singular, skip_singular,
      threshold_Lambda_ij);
  
  return new QuadraticPolynomialDOU(br_0, br_1, lengths, data);;
}

RCPP_EXPOSED_CLASS_NODECL(QuadraticPolynomialDOU::TraversalSpecificationType)
RCPP_EXPOSED_CLASS_NODECL(QuadraticPolynomialDOU::AlgorithmType)
  
RCPP_MODULE(PCMBaseCpp__QuadraticPolynomialDOU) {
  Rcpp::class_<QuadraticPolynomialDOU::TreeType::Tree> ( "PCMBaseCpp__QuadraticPolynomialDOU_Tree" )
  .property("num_nodes", &QuadraticPolynomialDOU::TreeType::Tree::num_nodes )
  .property("num_tips", &QuadraticPolynomialDOU::TreeType::Tree::num_tips )
  .method("FindNodeWithId", &QuadraticPolynomialDOU::TreeType::Tree::FindNodeWithId )
  .method("FindIdOfNode", &QuadraticPolynomialDOU::TreeType::Tree::FindIdOfNode )
  .method("FindIdOfParent", &QuadraticPolynomialDOU::TreeType::Tree::FindIdOfParent )
  .method("OrderNodes", &QuadraticPolynomialDOU::TreeType::Tree::OrderNodes )
  ;
  Rcpp::class_<QuadraticPolynomialDOU::TreeType>( "PCMBaseCpp__QuadraticPolynomialDOU_OrderedTree" )
    .derives<QuadraticPolynomialDOU::TreeType::Tree> ( "PCMBaseCpp__QuadraticPolynomialDOU_Tree" )
    .method("RangeIdPruneNode", &QuadraticPolynomialDOU::TreeType::RangeIdPruneNode )
    .method("RangeIdVisitNode", &QuadraticPolynomialDOU::TreeType::RangeIdVisitNode )
    .property("num_levels", &QuadraticPolynomialDOU::TreeType::num_levels )
    .property("ranges_id_visit", &QuadraticPolynomialDOU::TreeType::ranges_id_visit )
    .property("ranges_id_prune", &QuadraticPolynomialDOU::TreeType::ranges_id_prune )
  ;
  Rcpp::class_<QuadraticPolynomialDOU::AlgorithmType::ParentType>( "PCMBaseCpp__QuadraticPolynomialDOU_TraversalAlgorithm" )
    .property( "VersionOPENMP", &QuadraticPolynomialDOU::AlgorithmType::ParentType::VersionOPENMP )
    .property( "NumOmpThreads", &QuadraticPolynomialDOU::AlgorithmType::NumOmpThreads )
  ;
  Rcpp::class_<QuadraticPolynomialDOU::AlgorithmType> ( "PCMBaseCpp__QuadraticPolynomialDOU_ParallelPruning" )
    .derives<QuadraticPolynomialDOU::AlgorithmType::ParentType>( "PCMBaseCpp__QuadraticPolynomialDOU_TraversalAlgorithm" )
    .method( "ModeAutoStep", &QuadraticPolynomialDOU::AlgorithmType::ModeAutoStep )
    .property( "ModeAutoCurrent", &QuadraticPolynomialDOU::AlgorithmType::ModeAutoCurrent )
    .property( "IsTuning", &QuadraticPolynomialDOU::AlgorithmType::IsTuning )
    .property( "min_size_chunk_visit", &QuadraticPolynomialDOU::AlgorithmType::min_size_chunk_visit )
    .property( "min_size_chunk_prune", &QuadraticPolynomialDOU::AlgorithmType::min_size_chunk_prune )
    .property( "durations_tuning", &QuadraticPolynomialDOU::AlgorithmType::durations_tuning )
    .property( "fastest_step_tuning", &QuadraticPolynomialDOU::AlgorithmType::fastest_step_tuning )
  ;
  Rcpp::class_<QuadraticPolynomialDOU::TraversalSpecificationType::BaseType> ( "PCMBaseCpp__QuadraticPolynomial_PrunignSpec" )
    .field( "A", &QuadraticPolynomialDOU::TraversalSpecificationType::BaseType::A)
    .field( "b", &QuadraticPolynomialDOU::TraversalSpecificationType::BaseType::b )
    .field( "C", &QuadraticPolynomialDOU::TraversalSpecificationType::BaseType::C )
    .field( "d", &QuadraticPolynomialDOU::TraversalSpecificationType::BaseType::d )
    .field( "E", &QuadraticPolynomialDOU::TraversalSpecificationType::BaseType::E )
    .field( "f", &QuadraticPolynomialDOU::TraversalSpecificationType::BaseType::f )
    .field( "L", &QuadraticPolynomialDOU::TraversalSpecificationType::BaseType::L )
    .field( "m", &QuadraticPolynomialDOU::TraversalSpecificationType::BaseType::m )
    .field( "r", &QuadraticPolynomialDOU::TraversalSpecificationType::BaseType::r )
    .field( "V", &QuadraticPolynomialDOU::TraversalSpecificationType::BaseType::V )
    .field( "V_1", &QuadraticPolynomialDOU::TraversalSpecificationType::BaseType::V_1 )
  ;
  Rcpp::class_<QuadraticPolynomialDOU::TraversalSpecificationType> ( "PCMBaseCpp__QuadraticPolynomialDOU_PruningSpec" )
    .derives<QuadraticPolynomialDOU::TraversalSpecificationType::BaseType>( "PCMBaseCpp__QuadraticPolynomial_PrunignSpec" )
  ;
  Rcpp::class_<QuadraticPolynomialDOU>( "PCMBaseCpp__QuadraticPolynomialDOU" )
    .factory<arma::mat const&, Rcpp::List const&, Rcpp::List const&>(&CreateQuadraticPolynomialDOU)
    .method( "TraverseTree", &QuadraticPolynomialDOU::TraverseTree )
    .property( "tree", &QuadraticPolynomialDOU::tree )
    .property( "spec", &QuadraticPolynomialDOU::spec )
    .property( "algorithm", &QuadraticPolynomialDOU::algorithm )
  ;
}


QuadraticPolynomialMixedGaussian* CreateQuadraticPolynomialMixedGaussian(
    arma::mat const& X,
    Rcpp::List const& tree,
    Rcpp::List const& model,
    Rcpp::List const& metaInfo,
    std::vector<std::string> const& regimeModels) {
  
  double threshold_SV = static_cast<double>(metaInfo["PCMBase.Threshold.SV"]);
  double threshold_EV = static_cast<double>(metaInfo["PCMBase.Threshold.EV"]);
  double threshold_skip_singular = static_cast<double>(metaInfo["PCMBase.Threshold.Skip.Singular"]);
  double threshold_Lambda_ij = static_cast<double>(metaInfo["PCMBase.Threshold.Lambda_ij"]);
  
  bool skip_singular = static_cast<int>(metaInfo["PCMBase.Skip.Singular"]);
  
  arma::mat SE(Rcpp::as<arma::mat>(metaInfo["SE"]));
  
  Rcpp::List pcListInt = Rcpp::as<Rcpp::List>(metaInfo["pcListInt"]);
  std::vector<arma::uvec> Pc(Rcpp::as<arma::uword>(metaInfo["M"]));
  for(arma::uword i = 0; i < Pc.size(); ++i) {
    Pc[i] = Rcpp::as<arma::uvec>(pcListInt[i]);
  }
    
  arma::umat branches = tree["edge"];
  SPLITT::uvec br_0 = arma::conv_to<SPLITT::uvec>::from(branches.col(0));
  SPLITT::uvec br_1 = arma::conv_to<SPLITT::uvec>::from(branches.col(1));
  SPLITT::vec t = Rcpp::as<SPLITT::vec>(tree["edge.length"]);
  
  using namespace std;
  uint RModel = Rcpp::as<uint>(metaInfo["RModel"]);
  vector<arma::uword> regimes = Rcpp::as<vector<arma::uword> >(metaInfo["r"]);
  vector<arma::u8> jumps = Rcpp::as<vector<arma::u8> >(metaInfo["xi"]);
  
  if(regimes.size() != branches.n_rows) {
    ostringstream os;
    os<<"ERR:03841:PCMBaseCpp:Rcpp.cpp:CreateQuadraticPolynomialMixedGaussian:: The slot r in metaInfo has different length ("<<regimes.size()<<
      ") than the number of edges ("<<branches.n_rows<<").";
    throw logic_error(os.str());
  }
  
  if(jumps.size() != branches.n_rows) {
    ostringstream os;
    os<<"ERR:03842:PCMBaseCpp:Rcpp.cpp:CreateQuadraticPolynomialJOU:: The slot jumps in trees has different length ("<<jumps.size()<<
      ") than the number of edges ("<<branches.n_rows<<").";
    throw logic_error(os.str());
  }
  
  SPLITT::uint num_tips = Rcpp::as<Rcpp::CharacterVector>(tree["tip.label"]).size();
  
  SPLITT::uvec tip_names = SPLITT::Seq(static_cast<SPLITT::uint>(1), num_tips);
  
  vector<typename QuadraticPolynomialMixedGaussian::LengthType> lengths(branches.n_rows);
  
  for(arma::uword i = 0; i < branches.n_rows; ++i) {
    lengths[i].length_ = t[i];
    lengths[i].regime_ = regimes[i] - 1;
    lengths[i].jump_ = jumps[i];
  }
  
  
  if(threshold_SV <= 0) {
    ostringstream os;
    os<<"ERR:03843:PCMBaseCpp:Rcpp.cpp:CreateQuadraticPolynomialMixedGaussian:: The argument threshold_SV should be positive real number.";
    throw invalid_argument(os.str());
  }
  if(threshold_EV <= 0) {
    ostringstream os;
    os<<"ERR:03844:PCMBaseCpp:Rcpp.cpp:CreateQuadraticPolynomialMixedGaussian:: The argument threshold_EV should be positive real number.";
    throw invalid_argument(os.str());
  }
  if(threshold_Lambda_ij < 0) {
    ostringstream os;
    os<<"ERR:03845:PCMBaseCpp:Rcpp.cpp:CreateQuadraticPolynomialMixedGaussian:: The argument threshold_Lambda_ij should be non-negative double.";
    throw invalid_argument(os.str());
  }
  
  typename QuadraticPolynomialMixedGaussian::DataType data(
      tip_names, X, SE, Pc, 
      RModel, 
      regimeModels,
      threshold_SV, threshold_EV, threshold_skip_singular, skip_singular,
      threshold_Lambda_ij
      );
  
  return new QuadraticPolynomialMixedGaussian(br_0, br_1, lengths, data);
}

RCPP_EXPOSED_CLASS_NODECL(QuadraticPolynomialMixedGaussian::TraversalSpecificationType)
  RCPP_EXPOSED_CLASS_NODECL(QuadraticPolynomialMixedGaussian::AlgorithmType)
  
  RCPP_MODULE(PCMBaseCpp__QuadraticPolynomialMixedGaussian) {
    Rcpp::class_<QuadraticPolynomialMixedGaussian::TreeType::Tree> ( "PCMBaseCpp__QuadraticPolynomialMixedGaussian_Tree" )
    .property("num_nodes", &QuadraticPolynomialMixedGaussian::TreeType::Tree::num_nodes )
    .property("num_tips", &QuadraticPolynomialMixedGaussian::TreeType::Tree::num_tips )
    .method("FindNodeWithId", &QuadraticPolynomialMixedGaussian::TreeType::Tree::FindNodeWithId )
    .method("FindIdOfNode", &QuadraticPolynomialMixedGaussian::TreeType::Tree::FindIdOfNode )
    .method("FindIdOfParent", &QuadraticPolynomialMixedGaussian::TreeType::Tree::FindIdOfParent )
    .method("OrderNodes", &QuadraticPolynomialMixedGaussian::TreeType::Tree::OrderNodes )
    ;
    Rcpp::class_<QuadraticPolynomialMixedGaussian::TreeType>( "PCMBaseCpp__QuadraticPolynomialMixedGaussian_OrderedTree" )
      .derives<QuadraticPolynomialMixedGaussian::TreeType::Tree> ( "PCMBaseCpp__QuadraticPolynomialMixedGaussian_Tree" )
      .method("RangeIdPruneNode", &QuadraticPolynomialMixedGaussian::TreeType::RangeIdPruneNode )
      .method("RangeIdVisitNode", &QuadraticPolynomialMixedGaussian::TreeType::RangeIdVisitNode )
      .property("num_levels", &QuadraticPolynomialMixedGaussian::TreeType::num_levels )
      .property("ranges_id_visit", &QuadraticPolynomialMixedGaussian::TreeType::ranges_id_visit )
      .property("ranges_id_prune", &QuadraticPolynomialMixedGaussian::TreeType::ranges_id_prune )
    ;
    Rcpp::class_<QuadraticPolynomialMixedGaussian::AlgorithmType::ParentType>( "PCMBaseCpp__QuadraticPolynomialMixedGaussian_TraversalAlgorithm" )
      .property( "VersionOPENMP", &QuadraticPolynomialMixedGaussian::AlgorithmType::ParentType::VersionOPENMP )
      .property( "NumOmpThreads", &QuadraticPolynomialMixedGaussian::AlgorithmType::NumOmpThreads )
    ;
    Rcpp::class_<QuadraticPolynomialMixedGaussian::AlgorithmType> ( "PCMBaseCpp__QuadraticPolynomialMixedGaussian_ParallelPruning" )
      .derives<QuadraticPolynomialMixedGaussian::AlgorithmType::ParentType>( "PCMBaseCpp__QuadraticPolynomialMixedGaussian_TraversalAlgorithm" )
      .method( "ModeAutoStep", &QuadraticPolynomialMixedGaussian::AlgorithmType::ModeAutoStep )
      .property( "ModeAutoCurrent", &QuadraticPolynomialMixedGaussian::AlgorithmType::ModeAutoCurrent )
      .property( "IsTuning", &QuadraticPolynomialMixedGaussian::AlgorithmType::IsTuning )
      .property( "min_size_chunk_visit", &QuadraticPolynomialMixedGaussian::AlgorithmType::min_size_chunk_visit )
      .property( "min_size_chunk_prune", &QuadraticPolynomialMixedGaussian::AlgorithmType::min_size_chunk_prune )
      .property( "durations_tuning", &QuadraticPolynomialMixedGaussian::AlgorithmType::durations_tuning )
      .property( "fastest_step_tuning", &QuadraticPolynomialMixedGaussian::AlgorithmType::fastest_step_tuning )
    ;
    Rcpp::class_<QuadraticPolynomialMixedGaussian::TraversalSpecificationType::BaseType> ( "PCMBaseCpp__QuadraticPolynomial_PrunignSpec" )
      .field( "A", &QuadraticPolynomialMixedGaussian::TraversalSpecificationType::BaseType::A)
      .field( "b", &QuadraticPolynomialMixedGaussian::TraversalSpecificationType::BaseType::b )
      .field( "C", &QuadraticPolynomialMixedGaussian::TraversalSpecificationType::BaseType::C )
      .field( "d", &QuadraticPolynomialMixedGaussian::TraversalSpecificationType::BaseType::d )
      .field( "E", &QuadraticPolynomialMixedGaussian::TraversalSpecificationType::BaseType::E )
      .field( "f", &QuadraticPolynomialMixedGaussian::TraversalSpecificationType::BaseType::f )
      .field( "L", &QuadraticPolynomialMixedGaussian::TraversalSpecificationType::BaseType::L )
      .field( "m", &QuadraticPolynomialMixedGaussian::TraversalSpecificationType::BaseType::m )
      .field( "r", &QuadraticPolynomialMixedGaussian::TraversalSpecificationType::BaseType::r )
      .field( "V", &QuadraticPolynomialMixedGaussian::TraversalSpecificationType::BaseType::V )
      .field( "V_1", &QuadraticPolynomialMixedGaussian::TraversalSpecificationType::BaseType::V_1 )
    ;
    Rcpp::class_<QuadraticPolynomialMixedGaussian::TraversalSpecificationType> ( "PCMBaseCpp__QuadraticPolynomialMixedGaussian_PruningSpec" )
      .derives<QuadraticPolynomialMixedGaussian::TraversalSpecificationType::BaseType>( "PCMBaseCpp__QuadraticPolynomial_PrunignSpec" )
    ;
    Rcpp::class_<QuadraticPolynomialMixedGaussian>( "PCMBaseCpp__QuadraticPolynomialMixedGaussian" )
      .factory<arma::mat const&, Rcpp::List const&, Rcpp::List const&>(&CreateQuadraticPolynomialMixedGaussian)
      .method( "TraverseTree", &QuadraticPolynomialMixedGaussian::TraverseTree )
      .property( "tree", &QuadraticPolynomialMixedGaussian::tree )
      .property( "spec", &QuadraticPolynomialMixedGaussian::spec )
      .property( "algorithm", &QuadraticPolynomialMixedGaussian::algorithm )
    ;
  }
