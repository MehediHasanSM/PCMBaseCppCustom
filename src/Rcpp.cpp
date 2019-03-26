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

#include "QuadraticPolyWhite.h"
#include "QuadraticPolyBM.h"
#include "QuadraticPolyOU.h"
#include "QuadraticPolyJOU.h"
#include "QuadraticPolyDOU.h"
#include "QuadraticPolyMixedGaussian.h"


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

QuadraticPolyWhite* CreateQuadraticPolyWhite(
    arma::mat const&X, 
    Rcpp::List const& tree, 
    Rcpp::List const& model,
    Rcpp::List const& metaInfo) { 
  
  double threshold_SV = static_cast<double>(metaInfo["PCMBase.Threshold.SV"]);
  double threshold_EV = static_cast<double>(metaInfo["PCMBase.Threshold.EV"]);
  double threshold_skip_singular = static_cast<double>(metaInfo["PCMBase.Threshold.Skip.Singular"]);
  double threshold_Lambda_ij = static_cast<double>(metaInfo["PCMBase.Threshold.Lambda_ij"]);
  
  bool skip_singular = static_cast<int>(metaInfo["PCMBase.Skip.Singular"]);
  
  arma::cube VE(Rcpp::as<arma::cube>(metaInfo["VE"]));
    
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
    os<<"ERR:03801:PCMBaseCpp:Rcpp.cpp:CreateQuadraticPolyWhite:: The slot r in metaInfo has different length ("<<regimes.size()<<
      ") than the number of edges ("<<branches.n_rows<<").";
    throw logic_error(os.str());
  }
  
  SPLITT::uint num_tips = Rcpp::as<Rcpp::CharacterVector>(tree["tip.label"]).size();
  SPLITT::uvec tip_names = SPLITT::Seq(static_cast<SPLITT::uint>(1), num_tips);
  
  vector<typename QuadraticPolyWhite::LengthType> lengths(branches.n_rows);
  
  for(arma::uword i = 0; i < branches.n_rows; ++i) {
    lengths[i].length_ = t[i];
    lengths[i].regime_ = regimes[i] - 1;
  }
  
  
  if(threshold_SV <= 0) {
    ostringstream os;
    os<<"ERR:03802:PCMBaseCpp:Rcpp.cpp:CreateQuadraticPolyWhite:: The argument threshold_SV should be positive real number.";
    throw invalid_argument(os.str());
  }
  if(threshold_EV <= 0) {
    ostringstream os;
    os<<"ERR:03803:PCMBaseCpp:Rcpp.cpp:CreateQuadraticPolyMixedGaussian:: The argument threshold_EV should be positive real number.";
    throw invalid_argument(os.str());
  }
  typename QuadraticPolyWhite::DataType data(
      tip_names, X, VE, Pc, RModel, std::vector<std::string>(), 
      threshold_SV, threshold_EV, threshold_skip_singular, skip_singular,
      threshold_Lambda_ij);
  
  return new QuadraticPolyWhite(br_0, br_1, lengths, data);
}

  RCPP_EXPOSED_CLASS_NODECL(QuadraticPolyWhite::TreeType)
  RCPP_EXPOSED_CLASS_NODECL(QuadraticPolyWhite::TraversalSpecificationType)
  RCPP_EXPOSED_CLASS_NODECL(QuadraticPolyWhite::AlgorithmType)
  
  RCPP_MODULE(PCMBaseCpp__QuadraticPolyWhite) {
    Rcpp::class_<QuadraticPolyWhite::TreeType::Tree> ( "PCMBaseCpp__QuadraticPolyWhite_Tree" )
    .property("num_nodes", &QuadraticPolyWhite::TreeType::Tree::num_nodes )
    .property("num_tips", &QuadraticPolyWhite::TreeType::Tree::num_tips )
    .method("FindNodeWithId", &QuadraticPolyWhite::TreeType::Tree::FindNodeWithId )
    .method("FindIdOfNode", &QuadraticPolyWhite::TreeType::Tree::FindIdOfNode )
    .method("FindIdOfParent", &QuadraticPolyWhite::TreeType::Tree::FindIdOfParent )
    .method("OrderNodes", &QuadraticPolyWhite::TreeType::Tree::OrderNodes )
    ;
    Rcpp::class_<QuadraticPolyWhite::TreeType>( "PCMBaseCpp__QuadraticPolyWhite_OrderedTree" )
      .derives<QuadraticPolyWhite::TreeType::Tree> ( "PCMBaseCpp__QuadraticPolyWhite_Tree" )
      .method("RangeIdPruneNode", &QuadraticPolyWhite::TreeType::RangeIdPruneNode )
      .method("RangeIdVisitNode", &QuadraticPolyWhite::TreeType::RangeIdVisitNode )
      .property("num_levels", &QuadraticPolyWhite::TreeType::num_levels )
      .property("ranges_id_visit", &QuadraticPolyWhite::TreeType::ranges_id_visit )
      .property("ranges_id_prune", &QuadraticPolyWhite::TreeType::ranges_id_prune )
    ;
    Rcpp::class_<QuadraticPolyWhite::AlgorithmType::ParentType>( "PCMBaseCpp__QuadraticPolyWhite_TraversalAlgorithm" )
      .property( "VersionOPENMP", &QuadraticPolyWhite::AlgorithmType::ParentType::VersionOPENMP )
      .property( "NumOmpThreads", &QuadraticPolyWhite::AlgorithmType::NumOmpThreads )
    ;
    Rcpp::class_<QuadraticPolyWhite::AlgorithmType> ( "PCMBaseCpp__QuadraticPolyWhite_ParallelPruning" )
      .derives<QuadraticPolyWhite::AlgorithmType::ParentType>( "PCMBaseCpp__QuadraticPolyWhite_TraversalAlgorithm" )
      .method( "ModeAutoStep", &QuadraticPolyWhite::AlgorithmType::ModeAutoStep )
      .property( "ModeAutoCurrent", &QuadraticPolyWhite::AlgorithmType::ModeAutoCurrent )
      .property( "IsTuning", &QuadraticPolyWhite::AlgorithmType::IsTuning )
      .property( "min_size_chunk_visit", &QuadraticPolyWhite::AlgorithmType::min_size_chunk_visit )
      .property( "min_size_chunk_prune", &QuadraticPolyWhite::AlgorithmType::min_size_chunk_prune )
      .property( "durations_tuning", &QuadraticPolyWhite::AlgorithmType::durations_tuning )
      .property( "fastest_step_tuning", &QuadraticPolyWhite::AlgorithmType::fastest_step_tuning )
    ;
    Rcpp::class_<QuadraticPolyWhite::TraversalSpecificationType::BaseType> ( "PCMBaseCpp__QuadraticPoly_PrunignSpec" )
      .field( "A", &QuadraticPolyWhite::TraversalSpecificationType::BaseType::A)
      .field( "b", &QuadraticPolyWhite::TraversalSpecificationType::BaseType::b )
      .field( "C", &QuadraticPolyWhite::TraversalSpecificationType::BaseType::C )
      .field( "d", &QuadraticPolyWhite::TraversalSpecificationType::BaseType::d )
      .field( "E", &QuadraticPolyWhite::TraversalSpecificationType::BaseType::E )
      .field( "f", &QuadraticPolyWhite::TraversalSpecificationType::BaseType::f )
      .field( "L", &QuadraticPolyWhite::TraversalSpecificationType::BaseType::L )
      .field( "m", &QuadraticPolyWhite::TraversalSpecificationType::BaseType::m )
      .field( "r", &QuadraticPolyWhite::TraversalSpecificationType::BaseType::r )
      .field( "V", &QuadraticPolyWhite::TraversalSpecificationType::BaseType::V )
      .field( "V_1", &QuadraticPolyWhite::TraversalSpecificationType::BaseType::V_1 )
    ;
    Rcpp::class_<QuadraticPolyWhite::TraversalSpecificationType> ( "PCMBaseCpp__QuadraticPolyWhite_PruningSpec" )
      .derives<QuadraticPolyWhite::TraversalSpecificationType::BaseType>( "PCMBaseCpp__QuadraticPoly_PrunignSpec" )
    ;
    Rcpp::class_<QuadraticPolyWhite>( "PCMBaseCpp__QuadraticPolyWhite" )
      .factory<arma::mat const&, Rcpp::List const&, Rcpp::List const&>(&CreateQuadraticPolyWhite)
      .method( "TraverseTree", &QuadraticPolyWhite::TraverseTree )
      .property( "tree", &QuadraticPolyWhite::tree )
      .property( "spec", &QuadraticPolyWhite::spec )
      .property( "algorithm", &QuadraticPolyWhite::algorithm )
    ;
  }

QuadraticPolyBM* CreateQuadraticPolyBM(
    arma::mat const&X, 
    Rcpp::List const& tree, 
    Rcpp::List const& model,
    Rcpp::List const& metaInfo) { 
    
  double threshold_SV = static_cast<double>(metaInfo["PCMBase.Threshold.SV"]);
  double threshold_EV = static_cast<double>(metaInfo["PCMBase.Threshold.EV"]);
  double threshold_skip_singular = static_cast<double>(metaInfo["PCMBase.Threshold.Skip.Singular"]);
  double threshold_Lambda_ij = static_cast<double>(metaInfo["PCMBase.Threshold.Lambda_ij"]);
  
  bool skip_singular = static_cast<int>(metaInfo["PCMBase.Skip.Singular"]);
  
  arma::cube VE(Rcpp::as<arma::cube>(metaInfo["VE"]));
  
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
    os<<"ERR:03801:PCMBaseCpp:Rcpp.cpp:CreateQuadraticPolyBM:: The slot r in metaInfo has different length ("<<regimes.size()<<
      ") than the number of edges ("<<branches.n_rows<<").";
    throw logic_error(os.str());
  }
  
  SPLITT::uint num_tips = Rcpp::as<Rcpp::CharacterVector>(tree["tip.label"]).size();
  SPLITT::uvec tip_names = SPLITT::Seq(static_cast<SPLITT::uint>(1), num_tips);
  
  vector<typename QuadraticPolyBM::LengthType> lengths(branches.n_rows);
  
  for(arma::uword i = 0; i < branches.n_rows; ++i) {
    lengths[i].length_ = t[i];
    lengths[i].regime_ = regimes[i] - 1;
  }
  
  
  if(threshold_SV <= 0) {
    ostringstream os;
    os<<"ERR:03802:PCMBaseCpp:Rcpp.cpp:CreateQuadraticPolyBM:: The argument threshold_SV should be positive real number.";
    throw invalid_argument(os.str());
  }
  if(threshold_EV <= 0) {
    ostringstream os;
    os<<"ERR:03803:PCMBaseCpp:Rcpp.cpp:CreateQuadraticPolyMixedGaussian:: The argument threshold_EV should be positive real number.";
    throw invalid_argument(os.str());
  }
  typename QuadraticPolyBM::DataType data(
      tip_names, X, VE, Pc, RModel, std::vector<std::string>(), 
      threshold_SV, threshold_EV, threshold_skip_singular, skip_singular,
      threshold_Lambda_ij);
  
  return new QuadraticPolyBM(br_0, br_1, lengths, data);
}

//RCPP_EXPOSED_CLASS_NODECL(QuadraticPolyBM::TreeType)
RCPP_EXPOSED_CLASS_NODECL(QuadraticPolyBM::TraversalSpecificationType)
RCPP_EXPOSED_CLASS_NODECL(QuadraticPolyBM::AlgorithmType)
  
RCPP_MODULE(PCMBaseCpp__QuadraticPolyBM) {
  Rcpp::class_<QuadraticPolyBM::TreeType::Tree> ( "PCMBaseCpp__QuadraticPolyBM_Tree" )
  .property("num_nodes", &QuadraticPolyBM::TreeType::Tree::num_nodes )
  .property("num_tips", &QuadraticPolyBM::TreeType::Tree::num_tips )
  .method("FindNodeWithId", &QuadraticPolyBM::TreeType::Tree::FindNodeWithId )
  .method("FindIdOfNode", &QuadraticPolyBM::TreeType::Tree::FindIdOfNode )
  .method("FindIdOfParent", &QuadraticPolyBM::TreeType::Tree::FindIdOfParent )
  .method("OrderNodes", &QuadraticPolyBM::TreeType::Tree::OrderNodes )
  ;
  Rcpp::class_<QuadraticPolyBM::TreeType>( "PCMBaseCpp__QuadraticPolyBM_OrderedTree" )
    .derives<QuadraticPolyBM::TreeType::Tree> ( "PCMBaseCpp__QuadraticPolyBM_Tree" )
    .method("RangeIdPruneNode", &QuadraticPolyBM::TreeType::RangeIdPruneNode )
    .method("RangeIdVisitNode", &QuadraticPolyBM::TreeType::RangeIdVisitNode )
    .property("num_levels", &QuadraticPolyBM::TreeType::num_levels )
    .property("ranges_id_visit", &QuadraticPolyBM::TreeType::ranges_id_visit )
    .property("ranges_id_prune", &QuadraticPolyBM::TreeType::ranges_id_prune )
  ;
  Rcpp::class_<QuadraticPolyBM::AlgorithmType::ParentType>( "PCMBaseCpp__QuadraticPolyBM_TraversalAlgorithm" )
    .property( "VersionOPENMP", &QuadraticPolyBM::AlgorithmType::ParentType::VersionOPENMP )
    .property( "NumOmpThreads", &QuadraticPolyBM::AlgorithmType::NumOmpThreads )
  ;
  Rcpp::class_<QuadraticPolyBM::AlgorithmType> ( "PCMBaseCpp__QuadraticPolyBM_ParallelPruning" )
    .derives<QuadraticPolyBM::AlgorithmType::ParentType>( "PCMBaseCpp__QuadraticPolyBM_TraversalAlgorithm" )
    .method( "ModeAutoStep", &QuadraticPolyBM::AlgorithmType::ModeAutoStep )
    .property( "ModeAutoCurrent", &QuadraticPolyBM::AlgorithmType::ModeAutoCurrent )
    .property( "IsTuning", &QuadraticPolyBM::AlgorithmType::IsTuning )
    .property( "min_size_chunk_visit", &QuadraticPolyBM::AlgorithmType::min_size_chunk_visit )
    .property( "min_size_chunk_prune", &QuadraticPolyBM::AlgorithmType::min_size_chunk_prune )
    .property( "durations_tuning", &QuadraticPolyBM::AlgorithmType::durations_tuning )
    .property( "fastest_step_tuning", &QuadraticPolyBM::AlgorithmType::fastest_step_tuning )
  ;
  Rcpp::class_<QuadraticPolyBM::TraversalSpecificationType::BaseType> ( "PCMBaseCpp__QuadraticPoly_PrunignSpec" )
    .field( "A", &QuadraticPolyBM::TraversalSpecificationType::BaseType::A)
    .field( "b", &QuadraticPolyBM::TraversalSpecificationType::BaseType::b )
    .field( "C", &QuadraticPolyBM::TraversalSpecificationType::BaseType::C )
    .field( "d", &QuadraticPolyBM::TraversalSpecificationType::BaseType::d )
    .field( "E", &QuadraticPolyBM::TraversalSpecificationType::BaseType::E )
    .field( "f", &QuadraticPolyBM::TraversalSpecificationType::BaseType::f )
    .field( "L", &QuadraticPolyBM::TraversalSpecificationType::BaseType::L )
    .field( "m", &QuadraticPolyBM::TraversalSpecificationType::BaseType::m )
    .field( "r", &QuadraticPolyBM::TraversalSpecificationType::BaseType::r )
    .field( "V", &QuadraticPolyBM::TraversalSpecificationType::BaseType::V )
    .field( "V_1", &QuadraticPolyBM::TraversalSpecificationType::BaseType::V_1 )
  ;
  Rcpp::class_<QuadraticPolyBM::TraversalSpecificationType> ( "PCMBaseCpp__QuadraticPolyBM_PruningSpec" )
    .derives<QuadraticPolyBM::TraversalSpecificationType::BaseType>( "PCMBaseCpp__QuadraticPoly_PrunignSpec" )
  ;
  Rcpp::class_<QuadraticPolyBM>( "PCMBaseCpp__QuadraticPolyBM" )
    .factory<arma::mat const&, Rcpp::List const&, Rcpp::List const&>(&CreateQuadraticPolyBM)
    .method( "TraverseTree", &QuadraticPolyBM::TraverseTree )
    .property( "tree", &QuadraticPolyBM::tree )
    .property( "spec", &QuadraticPolyBM::spec )
    .property( "algorithm", &QuadraticPolyBM::algorithm )
  ;
}

QuadraticPolyOU* CreateQuadraticPolyOU(
    arma::mat const&X, 
    Rcpp::List const& tree, 
    Rcpp::List const& model,
    Rcpp::List const& metaInfo) {
  
  double threshold_SV = static_cast<double>(metaInfo["PCMBase.Threshold.SV"]);
  double threshold_EV = static_cast<double>(metaInfo["PCMBase.Threshold.EV"]);
  double threshold_skip_singular = static_cast<double>(metaInfo["PCMBase.Threshold.Skip.Singular"]);
  double threshold_Lambda_ij = static_cast<double>(metaInfo["PCMBase.Threshold.Lambda_ij"]);
  
  bool skip_singular = static_cast<int>(metaInfo["PCMBase.Skip.Singular"]);
  
  arma::cube VE(Rcpp::as<arma::cube>(metaInfo["VE"]));
  
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
    os<<"ERR:03811:PCMBaseCpp:Rcpp.cpp:CreateQuadraticPolyOU:: The slot r in metaInfo has different length ("<<regimes.size()<<
      ") than the number of edges ("<<branches.n_rows<<").";
    throw logic_error(os.str());
  }

  SPLITT::uint num_tips = Rcpp::as<Rcpp::CharacterVector>(tree["tip.label"]).size();
  SPLITT::uvec tip_names = SPLITT::Seq(static_cast<SPLITT::uint>(1), num_tips);
  
  vector<typename QuadraticPolyOU::LengthType> lengths(branches.n_rows);
  
  for(arma::uword i = 0; i < branches.n_rows; ++i) {
    lengths[i].length_ = t[i];
    lengths[i].regime_ = regimes[i] - 1;
  }
  
  if(threshold_SV <= 0) {
    ostringstream os;
    os<<"ERR:03812:PCMBaseCpp:Rcpp.cpp:CreateQuadraticPolyOU:: The argument threshold_SV should be positive real number.";
    throw invalid_argument(os.str());
  }
  if(threshold_EV <= 0) {
    ostringstream os;
    os<<"ERR:03814:PCMBaseCpp:Rcpp.cpp:CreateQuadraticPolyMixedGaussian:: The argument threshold_EV should be positive real number.";
    throw invalid_argument(os.str());
  }
  
  if(threshold_Lambda_ij < 0) {
    ostringstream os;
    os<<"ERR:03813:PCMBaseCpp:Rcpp.cpp:CreateQuadraticPolyOU:: The argument threshold_Lambda_ij should be non-negative double.";
    throw invalid_argument(os.str());
  }
  
  typename QuadraticPolyOU::DataType data(
      tip_names, X, VE, Pc, 
      RModel, std::vector<std::string>(), 
      threshold_SV, threshold_EV, threshold_skip_singular, skip_singular,
      threshold_Lambda_ij);
  
  return new QuadraticPolyOU(br_0, br_1, lengths, data);
}

RCPP_EXPOSED_CLASS_NODECL(QuadraticPolyOU::TraversalSpecificationType)
RCPP_EXPOSED_CLASS_NODECL(QuadraticPolyOU::AlgorithmType)
  
RCPP_MODULE(PCMBaseCpp__QuadraticPolyOU) {
  Rcpp::class_<QuadraticPolyOU::TreeType::Tree> ( "PCMBaseCpp__QuadraticPolyOU_Tree" )
  .property("num_nodes", &QuadraticPolyOU::TreeType::Tree::num_nodes )
  .property("num_tips", &QuadraticPolyOU::TreeType::Tree::num_tips )
  .method("FindNodeWithId", &QuadraticPolyOU::TreeType::Tree::FindNodeWithId )
  .method("FindIdOfNode", &QuadraticPolyOU::TreeType::Tree::FindIdOfNode )
  .method("FindIdOfParent", &QuadraticPolyOU::TreeType::Tree::FindIdOfParent )
  .method("OrderNodes", &QuadraticPolyOU::TreeType::Tree::OrderNodes )
  ;
  Rcpp::class_<QuadraticPolyOU::TreeType>( "PCMBaseCpp__QuadraticPolyOU_OrderedTree" )
    .derives<QuadraticPolyOU::TreeType::Tree> ( "PCMBaseCpp__QuadraticPolyOU_Tree" )
    .method("RangeIdPruneNode", &QuadraticPolyOU::TreeType::RangeIdPruneNode )
    .method("RangeIdVisitNode", &QuadraticPolyOU::TreeType::RangeIdVisitNode )
    .property("num_levels", &QuadraticPolyOU::TreeType::num_levels )
    .property("ranges_id_visit", &QuadraticPolyOU::TreeType::ranges_id_visit )
    .property("ranges_id_prune", &QuadraticPolyOU::TreeType::ranges_id_prune )
  ;
  Rcpp::class_<QuadraticPolyOU::AlgorithmType::ParentType>( "PCMBaseCpp__QuadraticPolyOU_TraversalAlgorithm" )
    .property( "VersionOPENMP", &QuadraticPolyOU::AlgorithmType::ParentType::VersionOPENMP )
    .property( "NumOmpThreads", &QuadraticPolyOU::AlgorithmType::NumOmpThreads )
  ;
  Rcpp::class_<QuadraticPolyOU::AlgorithmType> ( "PCMBaseCpp__QuadraticPolyOU_ParallelPruning" )
    .derives<QuadraticPolyOU::AlgorithmType::ParentType>( "PCMBaseCpp__QuadraticPolyOU_TraversalAlgorithm" )
    .method( "ModeAutoStep", &QuadraticPolyOU::AlgorithmType::ModeAutoStep )
    .property( "ModeAutoCurrent", &QuadraticPolyOU::AlgorithmType::ModeAutoCurrent )
    .property( "IsTuning", &QuadraticPolyOU::AlgorithmType::IsTuning )
    .property( "min_size_chunk_visit", &QuadraticPolyOU::AlgorithmType::min_size_chunk_visit )
    .property( "min_size_chunk_prune", &QuadraticPolyOU::AlgorithmType::min_size_chunk_prune )
    .property( "durations_tuning", &QuadraticPolyOU::AlgorithmType::durations_tuning )
    .property( "fastest_step_tuning", &QuadraticPolyOU::AlgorithmType::fastest_step_tuning )
  ;
  Rcpp::class_<QuadraticPolyOU::TraversalSpecificationType::BaseType> ( "PCMBaseCpp__QuadraticPoly_PrunignSpec" )
    .field( "A", &QuadraticPolyOU::TraversalSpecificationType::BaseType::A)
    .field( "b", &QuadraticPolyOU::TraversalSpecificationType::BaseType::b )
    .field( "C", &QuadraticPolyOU::TraversalSpecificationType::BaseType::C )
    .field( "d", &QuadraticPolyOU::TraversalSpecificationType::BaseType::d )
    .field( "E", &QuadraticPolyOU::TraversalSpecificationType::BaseType::E )
    .field( "f", &QuadraticPolyOU::TraversalSpecificationType::BaseType::f )
    .field( "L", &QuadraticPolyOU::TraversalSpecificationType::BaseType::L )
    .field( "m", &QuadraticPolyOU::TraversalSpecificationType::BaseType::m )
    .field( "r", &QuadraticPolyOU::TraversalSpecificationType::BaseType::r )
    .field( "V", &QuadraticPolyOU::TraversalSpecificationType::BaseType::V )
    .field( "V_1", &QuadraticPolyOU::TraversalSpecificationType::BaseType::V_1 )
  ;
  Rcpp::class_<QuadraticPolyOU::TraversalSpecificationType> ( "PCMBaseCpp__QuadraticPolyOU_PruningSpec" )
    .derives<QuadraticPolyOU::TraversalSpecificationType::BaseType>( "PCMBaseCpp__QuadraticPoly_PrunignSpec" )
  ;
  Rcpp::class_<QuadraticPolyOU>( "PCMBaseCpp__QuadraticPolyOU" )
    .factory<arma::mat const&, Rcpp::List const&, Rcpp::List const&>(&CreateQuadraticPolyOU)
    .method( "TraverseTree", &QuadraticPolyOU::TraverseTree )
    .property( "tree", &QuadraticPolyOU::tree )
    .property( "spec", &QuadraticPolyOU::spec )
    .property( "algorithm", &QuadraticPolyOU::algorithm )
  ;
}


QuadraticPolyJOU* CreateQuadraticPolyJOU(
    arma::mat const&X, 
    Rcpp::List const& tree, 
    Rcpp::List const& model,
    Rcpp::List const& metaInfo) {
  
  double threshold_SV = static_cast<double>(metaInfo["PCMBase.Threshold.SV"]);
  double threshold_EV = static_cast<double>(metaInfo["PCMBase.Threshold.EV"]);
  double threshold_skip_singular = static_cast<double>(metaInfo["PCMBase.Threshold.Skip.Singular"]);
  double threshold_Lambda_ij = static_cast<double>(metaInfo["PCMBase.Threshold.Lambda_ij"]);
  
  bool skip_singular = static_cast<int>(metaInfo["PCMBase.Skip.Singular"]);
 
  arma::cube VE(Rcpp::as<arma::cube>(metaInfo["VE"]));
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
    os<<"ERR:03821:PCMBaseCpp:Rcpp.cpp:CreateQuadraticPolyJOU:: The slot r in metaInfo has different length ("<<regimes.size()<<
      ") than the number of edges ("<<branches.n_rows<<").";
    throw logic_error(os.str());
  }
  
  if(jumps.size() != branches.n_rows) {
    ostringstream os;
    os<<"ERR:03822:PCMBaseCpp:Rcpp.cpp:CreateQuadraticPolyJOU:: The slot jumps in trees has different length ("<<jumps.size()<<
      ") than the number of edges ("<<branches.n_rows<<").";
    throw logic_error(os.str());
  }
  
  SPLITT::uint num_tips = Rcpp::as<Rcpp::CharacterVector>(tree["tip.label"]).size();
  SPLITT::uvec tip_names = SPLITT::Seq(static_cast<SPLITT::uint>(1), num_tips);
  
  vector<typename QuadraticPolyJOU::LengthType> lengths(branches.n_rows);
  
  for(arma::uword i = 0; i < branches.n_rows; ++i) {
    lengths[i].length_ = t[i];
    lengths[i].regime_ = regimes[i] - 1;
    lengths[i].jump_ = jumps[i];
  }
  
  if(threshold_SV <= 0) {
    ostringstream os;
    os<<"ERR:03823:PCMBaseCpp:Rcpp.cpp:CreateQuadraticPolyJOU:: The argument threshold_SV should be positive real number.";
    throw invalid_argument(os.str());
  }
  if(threshold_EV <= 0) {
    ostringstream os;
    os<<"ERR:03824:PCMBaseCpp:Rcpp.cpp:CreateQuadraticPolyMixedGaussian:: The argument threshold_EV should be positive real number.";
    throw invalid_argument(os.str());
  }
  if(threshold_Lambda_ij < 0) {
    ostringstream os;
    os<<"ERR:03825:PCMBaseCpp:Rcpp.cpp:CreateQuadraticPolyJOU:: The argument threshold_Lambda_ij should be non-negative double.";
    throw invalid_argument(os.str());
  }
  
  typename QuadraticPolyJOU::DataType data(
      tip_names, X, VE, Pc, 
      RModel, std::vector<std::string>(), 
      threshold_SV, threshold_EV, threshold_skip_singular, skip_singular,
      threshold_Lambda_ij);
  
  return new QuadraticPolyJOU(br_0, br_1, lengths, data);
}

RCPP_EXPOSED_CLASS_NODECL(QuadraticPolyJOU::TreeType)
RCPP_EXPOSED_CLASS_NODECL(QuadraticPolyJOU::TraversalSpecificationType)
RCPP_EXPOSED_CLASS_NODECL(QuadraticPolyJOU::AlgorithmType)
  
RCPP_MODULE(PCMBaseCpp__QuadraticPolyJOU) {
  Rcpp::class_<QuadraticPolyJOU::TreeType::Tree> ( "PCMBaseCpp__QuadraticPolyJOU_Tree" )
  .property("num_nodes", &QuadraticPolyJOU::TreeType::Tree::num_nodes )
  .property("num_tips", &QuadraticPolyJOU::TreeType::Tree::num_tips )
  .method("FindNodeWithId", &QuadraticPolyJOU::TreeType::Tree::FindNodeWithId )
  .method("FindIdOfNode", &QuadraticPolyJOU::TreeType::Tree::FindIdOfNode )
  .method("FindIdOfParent", &QuadraticPolyJOU::TreeType::Tree::FindIdOfParent )
  .method("OrderNodes", &QuadraticPolyJOU::TreeType::Tree::OrderNodes )
  ;
  Rcpp::class_<QuadraticPolyJOU::TreeType>( "PCMBaseCpp__QuadraticPolyJOU_OrderedTree" )
    .derives<QuadraticPolyJOU::TreeType::Tree> ( "PCMBaseCpp__QuadraticPolyJOU_Tree" )
    .method("RangeIdPruneNode", &QuadraticPolyJOU::TreeType::RangeIdPruneNode )
    .method("RangeIdVisitNode", &QuadraticPolyJOU::TreeType::RangeIdVisitNode )
    .property("num_levels", &QuadraticPolyJOU::TreeType::num_levels )
    .property("ranges_id_visit", &QuadraticPolyJOU::TreeType::ranges_id_visit )
    .property("ranges_id_prune", &QuadraticPolyJOU::TreeType::ranges_id_prune )
  ;
  Rcpp::class_<QuadraticPolyJOU::AlgorithmType::ParentType>( "PCMBaseCpp__QuadraticPolyJOU_TraversalAlgorithm" )
    .property( "VersionOPENMP", &QuadraticPolyJOU::AlgorithmType::ParentType::VersionOPENMP )
    .property( "NumOmpThreads", &QuadraticPolyJOU::AlgorithmType::NumOmpThreads )
  ;
  Rcpp::class_<QuadraticPolyJOU::AlgorithmType> ( "PCMBaseCpp__QuadraticPolyJOU_ParallelPruning" )
    .derives<QuadraticPolyJOU::AlgorithmType::ParentType>( "PCMBaseCpp__QuadraticPolyJOU_TraversalAlgorithm" )
    .method( "ModeAutoStep", &QuadraticPolyJOU::AlgorithmType::ModeAutoStep )
    .property( "ModeAutoCurrent", &QuadraticPolyJOU::AlgorithmType::ModeAutoCurrent )
    .property( "IsTuning", &QuadraticPolyJOU::AlgorithmType::IsTuning )
    .property( "min_size_chunk_visit", &QuadraticPolyJOU::AlgorithmType::min_size_chunk_visit )
    .property( "min_size_chunk_prune", &QuadraticPolyJOU::AlgorithmType::min_size_chunk_prune )
    .property( "durations_tuning", &QuadraticPolyJOU::AlgorithmType::durations_tuning )
    .property( "fastest_step_tuning", &QuadraticPolyJOU::AlgorithmType::fastest_step_tuning )
  ;
  Rcpp::class_<QuadraticPolyJOU::TraversalSpecificationType::BaseType> ( "PCMBaseCpp__QuadraticPoly_PrunignSpec" )
    .field( "A", &QuadraticPolyJOU::TraversalSpecificationType::BaseType::A)
    .field( "b", &QuadraticPolyJOU::TraversalSpecificationType::BaseType::b )
    .field( "C", &QuadraticPolyJOU::TraversalSpecificationType::BaseType::C )
    .field( "d", &QuadraticPolyJOU::TraversalSpecificationType::BaseType::d )
    .field( "E", &QuadraticPolyJOU::TraversalSpecificationType::BaseType::E )
    .field( "f", &QuadraticPolyJOU::TraversalSpecificationType::BaseType::f )
    .field( "L", &QuadraticPolyJOU::TraversalSpecificationType::BaseType::L )
    .field( "m", &QuadraticPolyJOU::TraversalSpecificationType::BaseType::m )
    .field( "r", &QuadraticPolyJOU::TraversalSpecificationType::BaseType::r )
    .field( "V", &QuadraticPolyJOU::TraversalSpecificationType::BaseType::V )
    .field( "V_1", &QuadraticPolyJOU::TraversalSpecificationType::BaseType::V_1 )
  ;
  Rcpp::class_<QuadraticPolyJOU::TraversalSpecificationType> ( "PCMBaseCpp__QuadraticPolyJOU_PruningSpec" )
    .derives<QuadraticPolyJOU::TraversalSpecificationType::BaseType>( "PCMBaseCpp__QuadraticPoly_PrunignSpec" )
  ;
  Rcpp::class_<QuadraticPolyJOU>( "PCMBaseCpp__QuadraticPolyJOU" )
    .factory<arma::mat const&, Rcpp::List const&, Rcpp::List const&>(&CreateQuadraticPolyJOU)
    .method( "TraverseTree", &QuadraticPolyJOU::TraverseTree )
    .property( "tree", &QuadraticPolyJOU::tree )
    .property( "spec", &QuadraticPolyJOU::spec )
    .property( "algorithm", &QuadraticPolyJOU::algorithm )
  ;
}

QuadraticPolyDOU* CreateQuadraticPolyDOU(
    arma::mat const& X,
    Rcpp::List const& tree,
    Rcpp::List const& model,
    Rcpp::List const& metaInfo) {
  
  double threshold_SV = static_cast<double>(metaInfo["PCMBase.Threshold.SV"]);
  double threshold_EV = static_cast<double>(metaInfo["PCMBase.Threshold.EV"]);
  double threshold_skip_singular = static_cast<double>(metaInfo["PCMBase.Threshold.Skip.Singular"]);
  double threshold_Lambda_ij = static_cast<double>(metaInfo["PCMBase.Threshold.Lambda_ij"]);
  
  bool skip_singular = static_cast<int>(metaInfo["PCMBase.Skip.Singular"]);

  arma::cube VE(Rcpp::as<arma::cube>(metaInfo["VE"]));
  
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
    os<<"ERR:03831:PCMBaseCpp:Rcpp.cpp:CreateQuadraticPolyDOU:: The slot r in metaInfo has different length ("<<regimes.size()<<
      ") than the number of edges ("<<branches.n_rows<<").";
    throw logic_error(os.str());
  }
  
  SPLITT::uint num_tips = Rcpp::as<Rcpp::CharacterVector>(tree["tip.label"]).size();
  SPLITT::uvec tip_names = SPLITT::Seq(static_cast<SPLITT::uint>(1), num_tips);
  
  vector<typename QuadraticPolyDOU::LengthType> lengths(branches.n_rows);
  
  for(arma::uword i = 0; i < branches.n_rows; ++i) {
    lengths[i].length_ = t[i];
    lengths[i].regime_ = regimes[i] - 1;
  }
  
  if(threshold_SV <= 0) {
    ostringstream os;
    os<<"ERR:03832:PCMBaseCpp:Rcpp.cpp:CreateQuadraticPolyDOU:: The argument threshold_SV should be positive real number.";
    throw invalid_argument(os.str());
  }
  if(threshold_EV <= 0) {
    ostringstream os;
    os<<"ERR:03834:PCMBaseCpp:Rcpp.cpp:CreateQuadraticPolyMixedGaussian:: The argument threshold_EV should be positive real number.";
    throw invalid_argument(os.str());
  }
  if(threshold_Lambda_ij < 0) {
    ostringstream os;
    os<<"ERR:03833:PCMBaseCpp:Rcpp.cpp:CreateQuadraticPolyDOU:: The argument threshold_Lambda_ij should be non-negative double.";
    throw invalid_argument(os.str());
  }
  
  typename QuadraticPolyDOU::DataType data(
      tip_names, X, VE, Pc, 
      RModel, std::vector<std::string>(), 
      threshold_SV, threshold_EV, threshold_skip_singular, skip_singular,
      threshold_Lambda_ij);
  
  return new QuadraticPolyDOU(br_0, br_1, lengths, data);;
}

RCPP_EXPOSED_CLASS_NODECL(QuadraticPolyDOU::TraversalSpecificationType)
RCPP_EXPOSED_CLASS_NODECL(QuadraticPolyDOU::AlgorithmType)
  
RCPP_MODULE(PCMBaseCpp__QuadraticPolyDOU) {
  Rcpp::class_<QuadraticPolyDOU::TreeType::Tree> ( "PCMBaseCpp__QuadraticPolyDOU_Tree" )
  .property("num_nodes", &QuadraticPolyDOU::TreeType::Tree::num_nodes )
  .property("num_tips", &QuadraticPolyDOU::TreeType::Tree::num_tips )
  .method("FindNodeWithId", &QuadraticPolyDOU::TreeType::Tree::FindNodeWithId )
  .method("FindIdOfNode", &QuadraticPolyDOU::TreeType::Tree::FindIdOfNode )
  .method("FindIdOfParent", &QuadraticPolyDOU::TreeType::Tree::FindIdOfParent )
  .method("OrderNodes", &QuadraticPolyDOU::TreeType::Tree::OrderNodes )
  ;
  Rcpp::class_<QuadraticPolyDOU::TreeType>( "PCMBaseCpp__QuadraticPolyDOU_OrderedTree" )
    .derives<QuadraticPolyDOU::TreeType::Tree> ( "PCMBaseCpp__QuadraticPolyDOU_Tree" )
    .method("RangeIdPruneNode", &QuadraticPolyDOU::TreeType::RangeIdPruneNode )
    .method("RangeIdVisitNode", &QuadraticPolyDOU::TreeType::RangeIdVisitNode )
    .property("num_levels", &QuadraticPolyDOU::TreeType::num_levels )
    .property("ranges_id_visit", &QuadraticPolyDOU::TreeType::ranges_id_visit )
    .property("ranges_id_prune", &QuadraticPolyDOU::TreeType::ranges_id_prune )
  ;
  Rcpp::class_<QuadraticPolyDOU::AlgorithmType::ParentType>( "PCMBaseCpp__QuadraticPolyDOU_TraversalAlgorithm" )
    .property( "VersionOPENMP", &QuadraticPolyDOU::AlgorithmType::ParentType::VersionOPENMP )
    .property( "NumOmpThreads", &QuadraticPolyDOU::AlgorithmType::NumOmpThreads )
  ;
  Rcpp::class_<QuadraticPolyDOU::AlgorithmType> ( "PCMBaseCpp__QuadraticPolyDOU_ParallelPruning" )
    .derives<QuadraticPolyDOU::AlgorithmType::ParentType>( "PCMBaseCpp__QuadraticPolyDOU_TraversalAlgorithm" )
    .method( "ModeAutoStep", &QuadraticPolyDOU::AlgorithmType::ModeAutoStep )
    .property( "ModeAutoCurrent", &QuadraticPolyDOU::AlgorithmType::ModeAutoCurrent )
    .property( "IsTuning", &QuadraticPolyDOU::AlgorithmType::IsTuning )
    .property( "min_size_chunk_visit", &QuadraticPolyDOU::AlgorithmType::min_size_chunk_visit )
    .property( "min_size_chunk_prune", &QuadraticPolyDOU::AlgorithmType::min_size_chunk_prune )
    .property( "durations_tuning", &QuadraticPolyDOU::AlgorithmType::durations_tuning )
    .property( "fastest_step_tuning", &QuadraticPolyDOU::AlgorithmType::fastest_step_tuning )
  ;
  Rcpp::class_<QuadraticPolyDOU::TraversalSpecificationType::BaseType> ( "PCMBaseCpp__QuadraticPoly_PrunignSpec" )
    .field( "A", &QuadraticPolyDOU::TraversalSpecificationType::BaseType::A)
    .field( "b", &QuadraticPolyDOU::TraversalSpecificationType::BaseType::b )
    .field( "C", &QuadraticPolyDOU::TraversalSpecificationType::BaseType::C )
    .field( "d", &QuadraticPolyDOU::TraversalSpecificationType::BaseType::d )
    .field( "E", &QuadraticPolyDOU::TraversalSpecificationType::BaseType::E )
    .field( "f", &QuadraticPolyDOU::TraversalSpecificationType::BaseType::f )
    .field( "L", &QuadraticPolyDOU::TraversalSpecificationType::BaseType::L )
    .field( "m", &QuadraticPolyDOU::TraversalSpecificationType::BaseType::m )
    .field( "r", &QuadraticPolyDOU::TraversalSpecificationType::BaseType::r )
    .field( "V", &QuadraticPolyDOU::TraversalSpecificationType::BaseType::V )
    .field( "V_1", &QuadraticPolyDOU::TraversalSpecificationType::BaseType::V_1 )
  ;
  Rcpp::class_<QuadraticPolyDOU::TraversalSpecificationType> ( "PCMBaseCpp__QuadraticPolyDOU_PruningSpec" )
    .derives<QuadraticPolyDOU::TraversalSpecificationType::BaseType>( "PCMBaseCpp__QuadraticPoly_PrunignSpec" )
  ;
  Rcpp::class_<QuadraticPolyDOU>( "PCMBaseCpp__QuadraticPolyDOU" )
    .factory<arma::mat const&, Rcpp::List const&, Rcpp::List const&>(&CreateQuadraticPolyDOU)
    .method( "TraverseTree", &QuadraticPolyDOU::TraverseTree )
    .property( "tree", &QuadraticPolyDOU::tree )
    .property( "spec", &QuadraticPolyDOU::spec )
    .property( "algorithm", &QuadraticPolyDOU::algorithm )
  ;
}


QuadraticPolyMixedGaussian* CreateQuadraticPolyMixedGaussian(
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
  
  arma::cube VE(Rcpp::as<arma::cube>(metaInfo["VE"]));
  
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
    os<<"ERR:03841:PCMBaseCpp:Rcpp.cpp:CreateQuadraticPolyMixedGaussian:: The slot r in metaInfo has different length ("<<regimes.size()<<
      ") than the number of edges ("<<branches.n_rows<<").";
    throw logic_error(os.str());
  }
  
  if(jumps.size() != branches.n_rows) {
    ostringstream os;
    os<<"ERR:03842:PCMBaseCpp:Rcpp.cpp:CreateQuadraticPolyJOU:: The slot jumps in trees has different length ("<<jumps.size()<<
      ") than the number of edges ("<<branches.n_rows<<").";
    throw logic_error(os.str());
  }
  
  SPLITT::uint num_tips = Rcpp::as<Rcpp::CharacterVector>(tree["tip.label"]).size();
  
  SPLITT::uvec tip_names = SPLITT::Seq(static_cast<SPLITT::uint>(1), num_tips);
  
  vector<typename QuadraticPolyMixedGaussian::LengthType> lengths(branches.n_rows);
  
  for(arma::uword i = 0; i < branches.n_rows; ++i) {
    lengths[i].length_ = t[i];
    lengths[i].regime_ = regimes[i] - 1;
    lengths[i].jump_ = jumps[i];
  }
  
  
  if(threshold_SV <= 0) {
    ostringstream os;
    os<<"ERR:03843:PCMBaseCpp:Rcpp.cpp:CreateQuadraticPolyMixedGaussian:: The argument threshold_SV should be positive real number.";
    throw invalid_argument(os.str());
  }
  if(threshold_EV <= 0) {
    ostringstream os;
    os<<"ERR:03844:PCMBaseCpp:Rcpp.cpp:CreateQuadraticPolyMixedGaussian:: The argument threshold_EV should be positive real number.";
    throw invalid_argument(os.str());
  }
  if(threshold_Lambda_ij < 0) {
    ostringstream os;
    os<<"ERR:03845:PCMBaseCpp:Rcpp.cpp:CreateQuadraticPolyMixedGaussian:: The argument threshold_Lambda_ij should be non-negative double.";
    throw invalid_argument(os.str());
  }
  
  typename QuadraticPolyMixedGaussian::DataType data(
      tip_names, X, VE, Pc, 
      RModel, 
      regimeModels,
      threshold_SV, threshold_EV, threshold_skip_singular, skip_singular,
      threshold_Lambda_ij
      );
  
  return new QuadraticPolyMixedGaussian(br_0, br_1, lengths, data);
}

RCPP_EXPOSED_CLASS_NODECL(QuadraticPolyMixedGaussian::TraversalSpecificationType)
  RCPP_EXPOSED_CLASS_NODECL(QuadraticPolyMixedGaussian::AlgorithmType)
  
  RCPP_MODULE(PCMBaseCpp__QuadraticPolyMixedGaussian) {
    Rcpp::class_<QuadraticPolyMixedGaussian::TreeType::Tree> ( "PCMBaseCpp__QuadraticPolyMixedGaussian_Tree" )
    .property("num_nodes", &QuadraticPolyMixedGaussian::TreeType::Tree::num_nodes )
    .property("num_tips", &QuadraticPolyMixedGaussian::TreeType::Tree::num_tips )
    .method("FindNodeWithId", &QuadraticPolyMixedGaussian::TreeType::Tree::FindNodeWithId )
    .method("FindIdOfNode", &QuadraticPolyMixedGaussian::TreeType::Tree::FindIdOfNode )
    .method("FindIdOfParent", &QuadraticPolyMixedGaussian::TreeType::Tree::FindIdOfParent )
    .method("OrderNodes", &QuadraticPolyMixedGaussian::TreeType::Tree::OrderNodes )
    ;
    Rcpp::class_<QuadraticPolyMixedGaussian::TreeType>( "PCMBaseCpp__QuadraticPolyMixedGaussian_OrderedTree" )
      .derives<QuadraticPolyMixedGaussian::TreeType::Tree> ( "PCMBaseCpp__QuadraticPolyMixedGaussian_Tree" )
      .method("RangeIdPruneNode", &QuadraticPolyMixedGaussian::TreeType::RangeIdPruneNode )
      .method("RangeIdVisitNode", &QuadraticPolyMixedGaussian::TreeType::RangeIdVisitNode )
      .property("num_levels", &QuadraticPolyMixedGaussian::TreeType::num_levels )
      .property("ranges_id_visit", &QuadraticPolyMixedGaussian::TreeType::ranges_id_visit )
      .property("ranges_id_prune", &QuadraticPolyMixedGaussian::TreeType::ranges_id_prune )
    ;
    Rcpp::class_<QuadraticPolyMixedGaussian::AlgorithmType::ParentType>( "PCMBaseCpp__QuadraticPolyMixedGaussian_TraversalAlgorithm" )
      .property( "VersionOPENMP", &QuadraticPolyMixedGaussian::AlgorithmType::ParentType::VersionOPENMP )
      .property( "NumOmpThreads", &QuadraticPolyMixedGaussian::AlgorithmType::NumOmpThreads )
    ;
    Rcpp::class_<QuadraticPolyMixedGaussian::AlgorithmType> ( "PCMBaseCpp__QuadraticPolyMixedGaussian_ParallelPruning" )
      .derives<QuadraticPolyMixedGaussian::AlgorithmType::ParentType>( "PCMBaseCpp__QuadraticPolyMixedGaussian_TraversalAlgorithm" )
      .method( "ModeAutoStep", &QuadraticPolyMixedGaussian::AlgorithmType::ModeAutoStep )
      .property( "ModeAutoCurrent", &QuadraticPolyMixedGaussian::AlgorithmType::ModeAutoCurrent )
      .property( "IsTuning", &QuadraticPolyMixedGaussian::AlgorithmType::IsTuning )
      .property( "min_size_chunk_visit", &QuadraticPolyMixedGaussian::AlgorithmType::min_size_chunk_visit )
      .property( "min_size_chunk_prune", &QuadraticPolyMixedGaussian::AlgorithmType::min_size_chunk_prune )
      .property( "durations_tuning", &QuadraticPolyMixedGaussian::AlgorithmType::durations_tuning )
      .property( "fastest_step_tuning", &QuadraticPolyMixedGaussian::AlgorithmType::fastest_step_tuning )
    ;
    Rcpp::class_<QuadraticPolyMixedGaussian::TraversalSpecificationType::BaseType> ( "PCMBaseCpp__QuadraticPoly_PrunignSpec" )
      .field( "A", &QuadraticPolyMixedGaussian::TraversalSpecificationType::BaseType::A)
      .field( "b", &QuadraticPolyMixedGaussian::TraversalSpecificationType::BaseType::b )
      .field( "C", &QuadraticPolyMixedGaussian::TraversalSpecificationType::BaseType::C )
      .field( "d", &QuadraticPolyMixedGaussian::TraversalSpecificationType::BaseType::d )
      .field( "E", &QuadraticPolyMixedGaussian::TraversalSpecificationType::BaseType::E )
      .field( "f", &QuadraticPolyMixedGaussian::TraversalSpecificationType::BaseType::f )
      .field( "L", &QuadraticPolyMixedGaussian::TraversalSpecificationType::BaseType::L )
      .field( "m", &QuadraticPolyMixedGaussian::TraversalSpecificationType::BaseType::m )
      .field( "r", &QuadraticPolyMixedGaussian::TraversalSpecificationType::BaseType::r )
      .field( "V", &QuadraticPolyMixedGaussian::TraversalSpecificationType::BaseType::V )
      .field( "V_1", &QuadraticPolyMixedGaussian::TraversalSpecificationType::BaseType::V_1 )
    ;
    Rcpp::class_<QuadraticPolyMixedGaussian::TraversalSpecificationType> ( "PCMBaseCpp__QuadraticPolyMixedGaussian_PruningSpec" )
      .derives<QuadraticPolyMixedGaussian::TraversalSpecificationType::BaseType>( "PCMBaseCpp__QuadraticPoly_PrunignSpec" )
    ;
    Rcpp::class_<QuadraticPolyMixedGaussian>( "PCMBaseCpp__QuadraticPolyMixedGaussian" )
      .factory<arma::mat const&, Rcpp::List const&, Rcpp::List const&>(&CreateQuadraticPolyMixedGaussian)
      .method( "TraverseTree", &QuadraticPolyMixedGaussian::TraverseTree )
      .property( "tree", &QuadraticPolyMixedGaussian::tree )
      .property( "spec", &QuadraticPolyMixedGaussian::spec )
      .property( "algorithm", &QuadraticPolyMixedGaussian::algorithm )
    ;
  }
