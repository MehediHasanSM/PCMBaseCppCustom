#ifndef QuadraticPoly_EB1D_H_
#define QuadraticPoly_EB1D_H_

#include "QuadraticPoly1D.h"
#include <armadillo>
#include <sstream>
#include <iostream>
#include <vector>

namespace PCMBaseCpp {

typedef SPLITT::OrderedTree<SPLITT::uint, LengthAndRegime> EB1DTreeType;

template<class TreeType, class DataType>
struct CondGaussianEB1D: public CondGaussianOmegaPhiV1D {
  
  TreeType const& ref_tree_;
  
  uint R_; 
  
  arma::vec X0;
  arma::vec Sigma;
  arma::vec Sigmae;
  arma::vec rho;
  
  // Pre-calculated node heights
  arma::vec height_;
  
  // Constructor 1: Takes R as an argument
  CondGaussianEB1D(TreeType const& ref_tree, DataType const& ref_data, uint R): ref_tree_(ref_tree) {
    this->R_ = R;
    
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
  CondGaussianEB1D(TreeType const& ref_tree, DataType const& ref_data): ref_tree_(ref_tree) {
    this->R_ = ref_data.R_;
    
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
    
    // Total params = X0(R) + rho(R) + Sigma_x(R) + Sigmae_x(R)
    uint npar = R_ * 4;
    if(par.size() - offset < npar) {
      std::ostringstream os;
      os << "QuadraticPolyEB1D.h:CondEB1D.SetParameter:: The length of the parameter vector minus offset (" << par.size() - offset <<
        ") should be at least of R*4, where R=" << R_ << " is the number of regimes.";
      throw std::logic_error(os.str());
    }
    
    X0 = vec(&par[offset], R_);
    rho = vec(&par[offset + R_], R_);
    Sigma = vec(&par[offset + 2*R_], R_);
    Sigmae = vec(&par[offset + 3*R_], R_);
    
    // By convention the parameters Sigma and Sigmae are passed as square roots
    for(uword r = 0; r < R_; r++) {
      Sigma(r) = Sigma(r) * Sigma(r);
      Sigmae(r) = Sigmae(r) * Sigmae(r);  
    }
    
    return npar;
  }
  
  void CalculateOmegaPhiV(uint i, arma::uword ri, arma::vec& omega, arma::vec& Phi, arma::vec& V) {
    using namespace arma;
    
    double ti = this->ref_tree_.LengthOfBranch(i).length_;
    double t_s = this->height_(this->ref_tree_.FindIdOfParent(i));
    
    omega(i) = 0.0;
    Phi(i) = 1.0;
    
    V(i) = ti * Sigma(ri) * exp(-rho(ri) * (ti + t_s) / 2.0);
    
    if(i < this->ref_tree_.num_tips()) {
      V(i) += Sigmae(ri);
    }
  }
};

class EB1D: public QuadraticPoly1D<EB1DTreeType> {
public:
  typedef EB1DTreeType TreeType;
  typedef QuadraticPoly1D<TreeType> BaseType;
  typedef EB1D MyType;
  typedef arma::vec StateType;
  typedef NumericTraitData1D<TreeType::NodeType> DataType;
  typedef std::vector<double> ParameterType;
  typedef SPLITT::PostOrderTraversal<MyType> AlgorithmType;
  
  CondGaussianEB1D<TreeType, DataType> cond_dist_;
  
  EB1D(TreeType const& tree, DataType const& input_data):
    BaseType(tree, input_data), cond_dist_(tree, input_data) {
    BaseType::ptr_cond_dist_.push_back(&cond_dist_);
  }
  
  void SetParameter(ParameterType const& par) {
    cond_dist_.SetParameter(par, 0);
  }
};

typedef TraversalTaskWrapper<EB1D> QuadraticPolyEB1D;

} // namespace PCMBaseCpp

#endif // QuadraticPoly_EB1D_H_