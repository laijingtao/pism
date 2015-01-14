// Copyright (C) 2012, 2014, 2015  David Maxwell and Constantine Khroulev
//
// This file is part of PISM.
//
// PISM is free software; you can redistribute it and/or modify it under the
// terms of the GNU General Public License as published by the Free Software
// Foundation; either version 3 of the License, or (at your option) any later
// version.
//
// PISM is distributed in the hope that it will be useful, but WITHOUT ANY
// WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
// FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
// details.
//
// You should have received a copy of the GNU General Public License
// along with PISM; if not, write to the Free Software
// Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA

#include "IP_H1NormFunctional.hh"
#include "error_handling.hh"
#include "IceGrid.hh"

namespace pism {

PetscErrorCode IP_H1NormFunctional2S::valueAt(IceModelVec2S &x, double *OUTPUT) {

  // The value of the objective
  double value = 0;

  double x_e[FEQuadrature::Nk];
  double x_q[FEQuadrature::Nq], dxdx_q[FEQuadrature::Nq], dxdy_q[FEQuadrature::Nq];
  IceModelVec::AccessList list(x);

  // Jacobian times weights for quadrature.
  const double* JxW = m_quadrature.getWeightedJacobian();

  DirichletData_Scalar dirichletBC;
  dirichletBC.init(m_dirichletIndices, NULL);

  // Loop through all LOCAL elements.
  int xs = m_element_index.lxs, xm = m_element_index.lxm,
           ys = m_element_index.lys, ym = m_element_index.lym;
  for (int i=xs; i<xs+xm; i++) {
    for (int j=ys; j<ys+ym; j++) {
      m_dofmap.reset(i, j, m_grid);

      // Obtain values of x at the quadrature points for the element.
      m_dofmap.extractLocalDOFs(x, x_e);
      if (dirichletBC) {
        dirichletBC.update_homogeneous(m_dofmap, x_e);
      }
      m_quadrature.computeTrialFunctionValues(x_e, x_q, dxdx_q, dxdy_q);

      for (unsigned int q=0; q<FEQuadrature::Nq; q++) {
        value += JxW[q]*(m_cL2*x_q[q]*x_q[q]+ m_cH1*(dxdx_q[q]*dxdx_q[q]+dxdy_q[q]*dxdy_q[q]));
      } // q
    } // j
  } // i

  GlobalSum(m_grid.com, &value, OUTPUT, 1);

  dirichletBC.finish();


  return 0;
}

PetscErrorCode IP_H1NormFunctional2S::dot(IceModelVec2S &a, IceModelVec2S &b, double *OUTPUT) {

  // The value of the objective
  double value = 0;

  double a_e[FEQuadrature::Nk];
  double a_q[FEQuadrature::Nq], dadx_q[FEQuadrature::Nq], dady_q[FEQuadrature::Nq];

  double b_e[FEQuadrature::Nk];
  double b_q[FEQuadrature::Nq], dbdx_q[FEQuadrature::Nq], dbdy_q[FEQuadrature::Nq];

  IceModelVec::AccessList list(a);
  list.add(b);

  // Jacobian times weights for quadrature.
  const double* JxW = m_quadrature.getWeightedJacobian();

  DirichletData_Scalar dirichletBC;
  dirichletBC.init(m_dirichletIndices, NULL);

  // Loop through all LOCAL elements.
  int xs = m_element_index.lxs, xm = m_element_index.lxm,
           ys = m_element_index.lys, ym = m_element_index.lym;
  for (int i=xs; i<xs+xm; i++) {
    for (int j=ys; j<ys+ym; j++) {
      m_dofmap.reset(i, j, m_grid);

      // Obtain values of x at the quadrature points for the element.
      m_dofmap.extractLocalDOFs(a, a_e);
      if (dirichletBC) {
        dirichletBC.update_homogeneous(m_dofmap, a_e);
      }
      m_quadrature.computeTrialFunctionValues(a_e, a_q, dadx_q, dady_q);

      m_dofmap.extractLocalDOFs(b, b_e);
      if (dirichletBC) {
        dirichletBC.update_homogeneous(m_dofmap, b_e);
      }
      m_quadrature.computeTrialFunctionValues(b_e, b_q, dbdx_q, dbdy_q);

      for (unsigned int q=0; q<FEQuadrature::Nq; q++) {
        value += JxW[q]*(m_cL2*a_q[q]*b_q[q]+ m_cH1*(dadx_q[q]*dbdx_q[q]+dady_q[q]*dbdy_q[q]));
      } // q
    } // j
  } // i

  GlobalSum(m_grid.com, &value, OUTPUT, 1);

  dirichletBC.finish();

  return 0;
}


PetscErrorCode IP_H1NormFunctional2S::gradientAt(IceModelVec2S &x, IceModelVec2S &gradient) {

  // Clear the gradient before doing anything with it!
  gradient.set(0);

  double x_e[FEQuadrature::Nk];
  double x_q[FEQuadrature::Nq], dxdx_q[FEQuadrature::Nq], dxdy_q[FEQuadrature::Nq];

  double gradient_e[FEQuadrature::Nk];

  IceModelVec::AccessList list(x);
  list.add(gradient);

  // An Nq by Nk array of test function values.
  const FEFunctionGerm (*test)[FEQuadrature::Nk] = m_quadrature.testFunctionValues();

  // Jacobian times weights for quadrature.
  const double* JxW = m_quadrature.getWeightedJacobian();

  DirichletData_Scalar dirichletBC;
  dirichletBC.init(m_dirichletIndices, NULL);

  // Loop through all local and ghosted elements.
  int xs = m_element_index.xs, xm = m_element_index.xm,
           ys = m_element_index.ys, ym = m_element_index.ym;
  for (int i=xs; i<xs+xm; i++) {
    for (int j=ys; j<ys+ym; j++) {

      // Reset the DOF map for this element.
      m_dofmap.reset(i, j, m_grid);

      // Obtain values of x at the quadrature points for the element.
      m_dofmap.extractLocalDOFs(i, j, x, x_e);
      if (dirichletBC) {
        dirichletBC.constrain(m_dofmap);
        dirichletBC.update_homogeneous(m_dofmap, x_e);
      }
      m_quadrature.computeTrialFunctionValues(x_e, x_q, dxdx_q, dxdy_q);

      // Zero out the element-local residual in prep for updating it.
      for (unsigned int k=0; k<FEQuadrature::Nk; k++) {
        gradient_e[k] = 0;
      }

      for (unsigned int q=0; q<FEQuadrature::Nq; q++) {
        const double &x_qq=x_q[q];
        const double &dxdx_qq=dxdx_q[q], &dxdy_qq=dxdy_q[q];
        for (unsigned int k=0; k<FEQuadrature::Nk; k++) {
          gradient_e[k] += 2*JxW[q]*(m_cL2*x_qq*test[q][k].val +
            m_cH1*(dxdx_qq*test[q][k].dx + dxdy_qq*test[q][k].dy));
        } // k
      } // q
      m_dofmap.addLocalResidualBlock(gradient_e, gradient);
    } // j
  } // i

  dirichletBC.finish();

  return 0;
}

PetscErrorCode IP_H1NormFunctional2S::assemble_form(Mat form) {
  PetscErrorCode   ierr;

  // Zero out the Jacobian in preparation for updating it.
  MatZeroEntries(form);

  // Jacobian times weights for quadrature.
  const double* JxW = m_quadrature.getWeightedJacobian();

  DirichletData_Scalar zeroLocs;
  zeroLocs.init(m_dirichletIndices, NULL);

  // Values of the finite element test functions at the quadrature points.
  // This is an Nq by Nk array of function germs (Nq=#of quad pts, Nk=#of test functions).
  const FEFunctionGerm (*test)[FEQuadrature::Nk] = m_quadrature.testFunctionValues();

  // Loop through all the elements.
  int xs = m_element_index.xs,
    xm   = m_element_index.xm,
    ys   = m_element_index.ys,
    ym   = m_element_index.ym;
  for (int i=xs; i<xs+xm; i++) {
    for (int j=ys; j<ys+ym; j++) {
      // Element-local Jacobian matrix (there are FEQuadrature::Nk vector valued degrees
      // of freedom per elment, for a total of (2*FEQuadrature::Nk)*(2*FEQuadrature::Nk) = 16
      // entries in the local Jacobian.
      double K[FEQuadrature::Nk][FEQuadrature::Nk];

      // Initialize the map from global to local degrees of freedom for this element.
      m_dofmap.reset(i, j, m_grid);

      // Don't update rows/cols where we project to zero.
      if (zeroLocs) {
        zeroLocs.constrain(m_dofmap);
      }

      // Build the element-local Jacobian.
      ierr = PetscMemzero(K, sizeof(K));
      PISM_PETSC_CHK(ierr, "PetscMemzero");
      for (unsigned int q=0; q<FEQuadrature::Nq; q++) {
        for (unsigned int k = 0; k < FEQuadrature::Nk; k++) {   // Test functions
          for (unsigned int l = 0; l < FEQuadrature::Nk; l++) { // Trial functions
            const FEFunctionGerm &test_qk=test[q][k];
            const FEFunctionGerm &test_ql=test[q][l];
            K[k][l] += JxW[q]*(m_cL2*test_qk.val*test_ql.val +
                               m_cH1*(test_qk.dx*test_ql.dx +
                                      test_qk.dy*test_ql.dy));
          } // l
        } // k
      } // q
      m_dofmap.addLocalJacobianBlock(&K[0][0], form);
    } // j
  } // i

  if (zeroLocs) {
    zeroLocs.fix_jacobian(form);
  }
  zeroLocs.finish();

  MatAssemblyBegin(form, MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(form, MAT_FINAL_ASSEMBLY);

  return 0;
}

} // end of namespace pism
