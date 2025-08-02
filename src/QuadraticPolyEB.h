#ifndef QuadraticPoly_EB_H_
#define QuadraticPoly_EB_H_

#include "QuadraticPoly.h"
#include <armadillo>
#include <sstream>
#include <iostream>
#include <vector>

namespace PCMBaseCpp {

typedef SPLITT::OrderedTree<SPLITT::uint, LengthAndRegime> EBTreeType;

template<class TreeType, class DataType>
struct CondGaussianEB: public CondGaussianOmegaPhiV {
  
  TreeType const& ref_tree_;
  
  uint k_;
  uint R_; 
  bool transpose_Sigma_x = false; 
  
  arma::mat X0;
  arma::cube Sigma;
  arma::cube Sigmae;
  arma::vec rho;
  arma::mat I;
  
  // variable to store pre-calculated node heights
  arma::vec height_;
  
  // Constructor 1: Takes R as an argument
  CondGaussianEB(TreeType const& ref_tree, DataType const& ref_data, uint R): ref_tree_(ref_tree) {
    this->k_ = ref_data.k_;
    this->R_ = R;
    this->I = arma::eye(k_, k_);
    this->transpose_Sigma_x = ref_data.transpose_Sigma_x;
    
    // Pre-calculate node heights
    height_.set_size(ref_tree_.num_nodes());
    std::vector<uint> q;
    q.push_back(ref_tree_.num_nodes() - 1); // Start with the root
    height_(ref_tree_.num_nodes() - 1) = 0.0;
    
    int head = 0;
    while(head < q.size()) {
      uint p_id = q[head];
      head++;
      for(uint c_id : ref_tree_.FindChildren(p_id)) {
        height_(c_id) = height_(p_id) + ref_tree_.LengthOfBranch(c_id).length_;
        q.push_back(c_id);
      }
    }
  }
  
  // Constructor 2: Infers R from the data object
  CondGaussianEB(TreeType const& ref_tree, DataType const& ref_data): ref_tree_(ref_tree) {
    this->k_ = ref_data.k_;
    this->R_ = ref_data.R_;
    this->I = arma::eye(k_, k_);
    this->transpose_Sigma_x = ref_data.transpose_Sigma_x;
    
    // Pre-calculate node heights
    height_.set_size(ref_tree_.num_nodes());
    std::vector<uint> q;
    q.push_back(ref_tree_.num_nodes() - 1); // Start with the root
    height_(ref_tree_.num_nodes() - 1) = 0.0;
    
    int head = 0;
    while(head < q.size()) {
      uint p_id = q[head];
      head++;
      for(uint c_id : ref_tree_.FindChildren(p_id)) {
        height_(c_id) = height_(p_id) + ref_tree_.LengthOfBranch(c_id).length_;
        q.push_back(c_id);
      }
    }
  }
  
  arma::uword SetParameter(std::vector<double> const& par, arma::uword offset) {
    using namespace arma;
    
    uint npar = R_ * (k_ + 1 + 2*k_*k_);
    if(par.size() - offset < npar) {
      std::ostringstream os;
      os << "QuadraticPolyEB.h:CondEB.SetParameter:: The length of the parameter vector minus offset (" << par.size() - offset <<
        ") should be at least of R*(k+1+2*k^2), where k=" << k_ << " is the number of traits and " <<
          " R=" << R_ << " is the number of regimes.";
      throw std::logic_error(os.str());
    }
    
    X0 = mat(&par[offset], k_, R_);
    rho = vec(&par[offset + k_*R_], R_);
    Sigma = cube(&par[offset + k_*R_ + R_], k_, k_, R_);
    Sigmae = cube(&par[offset + k_*R_ + R_ + k_*k_*R_], k_, k_, R_);
    
    if(transpose_Sigma_x) {
      for(uword r = 0; r < R_; r++) {
        Sigma.slice(r) = Sigma.slice(r).t() * Sigma.slice(r);
        Sigmae.slice(r) = Sigmae.slice(r).t() * Sigmae.slice(r);  
      }
    } else {
      for(uword r = 0; r < R_; r++) {
        Sigma.slice(r) = Sigma.slice(r) * Sigma.slice(r).t();
        Sigmae.slice(r) = Sigmae.slice(r) * Sigmae.slice(r).t();  
      }
    }
    return npar;
  }
  
  void CalculateOmegaPhiV(uint i, arma::uword ri, arma::mat& omega, arma::cube& Phi, arma::cube& V) {
    using namespace arma;
    
    // branch length (equivalent to 't' in R)
    double ti = this->ref_tree_.LengthOfBranch(i).length_;
    
    // time from root to the start of the current branch (equivalent to 't_s' in R)
    double t_s = this->height_(this->ref_tree_.FindIdOfParent(i));
    
    omega.col(i).fill(0);
    Phi.slice(i) = I;
    
    V.slice(i) = ti * Sigma.slice(ri) * exp(-rho(ri) * (ti + t_s) / 2.0); 
    
    // Add non-heritable variance for tips
    if(i < this->ref_tree_.num_tips()) {
      V.slice(i) += Sigmae.slice(ri);
    }
  }
};

class EB: public QuadraticPoly<EBTreeType> {
public:
  typedef EBTreeType TreeType;
  typedef QuadraticPoly<TreeType> BaseType;
  typedef EB MyType;
  typedef arma::vec StateType;
  typedef NumericTraitData<TreeType::NodeType> DataType;
  typedef std::vector<double> ParameterType;
  typedef SPLITT::PostOrderTraversal<MyType> AlgorithmType;
  
  CondGaussianEB<TreeType, DataType> cond_dist_;
  
  EB(TreeType const& tree, DataType const& input_data):
    BaseType(tree, input_data), cond_dist_(tree, input_data) {
    BaseType::ptr_cond_dist_.push_back(&cond_dist_);
  }
  
  void SetParameter(ParameterType const& par) {
    cond_dist_.SetParameter(par, 0);
  }
};

typedef TraversalTaskWrapper<EB> QuadraticPolyEB;

} // namespace PCMBaseCpp

#endif // QuadraticPoly_EB_H_