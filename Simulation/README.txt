#===============================================================================================
Description: This folder contains all data, scripts, and write ups associated with the Simulations of the MRGN, MRPC, and GMAC for trios with
all types of confounding variables (confounder, intermediate, common child). see /sim_write_up_files/ folder and check the file with latest date for details
regarding the specifics of the simulations.

Other Notes:

*the folder /adapted_GMAC_func/ contains slightly modified source code for GMAC that allows for timing of the confounder selection and inference steps --> which is needed for the 
 computation time plots in the MRGN manuscript see /MRGN_extra/Manuscript/ for details

*all simulated trios are in the /sim_set_3_data/ folder which is organized as follows:

                                 sim_set_3_data
                                       |
                                       |
          ---------------------------------------------------------------------------------------------------------
          |                                 |                                 |                                   |
          |                                 |                                 |                                   |
     diagnostics                int_and_child_filtered_data             many_conf_data                support_data_for_write_ups
  (comp time data)              (simulations with < 15confs)       (simulations with 15-50 confs)    (other tables/data from write ups)
          |                                 |                                 |                                      
          |                                 |                                 |          
         \ /                               \ /                               \ /
          *                                 *                                 *
          ---------------------------------------------------------------------
                                            |
                                            |
                                           \ /
                                            *
                               alpha_01_selected_confs_results
                        -----------------------------------------------
                       | results for MRGN when using alpha<01 cuttoff  | 
                       | to select the confounders - in place of usual |
                       | FDR <10%                                      |
                        -----------------------------------------------






