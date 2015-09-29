#ReadMe file of sol_package directory#

#Structure of sol_package#
sol_package/talk       /k.osaki_msu_seminar_9_25_2015.pptx
					   /figures/basic_plot/
					           /convergence_of_conserved_quantities/
					   		   /large_ptheta/
					  		   /Theory_order_conf/
					 		   /Mag_field_of_several_sol_length/
							   /Ptheta_0_sol4cm/
							   /Ptheta_0_sol30cm/
							   /Ptheta_finit_sol4cm/
							   /Ptheta_finit_sol30cm/					   
		   /orbit_solve/sol_scipy_sample.sh
		   			   /sol_scipy_sample.py
					   /energy_cal_sample.nb
					   /sample_scipy_nonlinear_sol_basic plot.nb

#Explanation of each directory and file#
There are two directories in sol_package "talk" and "orbit_solve".
In the talk directory, there are the powerpoint_file of osaki's seminar and figure directory.
In the figure directory, there are nine directories.
  1.(*)basic_plot : We can make plots of conserved quantities and orbit plots by sample_scipy_nonlinear_sol_basic plot.nb.
  2.(*)convergence_of_conserved_quantities : We can make plots of conserved quantities with several step size by sample_scipy_nonlinear_sol_nstep.nb.
  3.large_ptheta : We can calculate upper limit of P_theta with several Brho by Confirm_large_Ptheta_sample.nb.
  4.Theory_order_conf : We can see the magnitudes of theoritical radial kick of each order by Thoeretical_cal_order1.nb.
  5.Mag_field_of_several_sol_length : We can make plots of magnetic fields and derivation of Bz by B_field_cal_several_aspect_ratio.nb and B_field_cal_sample.nb.
  6.(*)Ptheta_0_sol4cm : We can see the difference of radial kick between theory and simulation in the case of Ptheta=0 and sol=4cm by sample_scipy_nonlinear_sol4cm_pth0.nb.
  7.(*)Ptheta_0_sol30cm : We can see the difference of radial kick between theory and simulation in the case of Ptheta=0 and sol=30cm by sample_scipy_nonlinear_sol30cm_pth0.nb.
  8.(*)Ptheta_finit_sol4cm : We can see the difference of radial kick between theory and simulation in the case of Ptheta!=0 and sol=4cm by sample_scipy_nonlinear_sol4cm_pth_finit.nb.
  9.(*)Ptheta_finit_sol30cm : We can see the difference of radial kick between theory and simulation in the case of Ptheta!=0 and sol=30cm by sample_scipy_nonlinear_sol30cm_pth_finit.nb.
In the case of (*) directories, the proper data files are needed.

In the orbit_solve directory, there are sample scipy programs of a single particle simulation. We can conduct simulation by "./sol_scipy_sample.sh".
From the energy_cal_sample.nb, we can calculate beam energies and we can use them in sol_scipy_sample.sh.
After conduct "./sol_scipy_sample.sh", we obtain orbit data and can make plots by sample_scipy_nonlinear_sol_basic plot.nb.