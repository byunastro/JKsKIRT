import numpy as np

import os
import subprocess

from . import make_INPUT
from . import make_SKI

#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

def execute(repo,
            repo_output,
            ski_name,
            skirt_dir,
            N_thread):
    
    os.chdir(repo+repo_output)
    process = subprocess.Popen("%s/skirt -t %s ../%s"%(skirt_dir,N_thread,ski_name),
                               shell=True,
                               universal_newlines=True,
                               stdout=subprocess.PIPE)
    for line in process.stdout:
        print(line, end='')

    return 


def delete_indat(repo,
                 opart_name,
                 ypart_name,
                 gas_name):

    return 



def make_INSKI(
               boxtokpc,
               x_s,y_s,z_s,
               vx_s,vy_s,vz_s,
               m_s,m0_s,
               age_s,metal_s,
               x_c,y_c,z_c,lvl_c,
               vx_c,vy_c,vz_c,
               m_c,
               T_c,metal_c,
               param,
               pos_ctr,
               vel_ctr,    
                ):
    
    os.makedirs(param.repo+param.repo_output,exist_ok=True)
    
    make_INPUT.make_input(
                           boxtokpc=boxtokpc,
                           x_s=x_s,
                           y_s=y_s,
                           z_s=z_s,
                           vx_s=vx_s,
                           vy_s=vy_s,
                           vz_s=vz_s,
                           m_s=m_s,
                           m0_s=m0_s,
                           age_s=age_s,
                           metal_s=metal_s,
                           x_c=x_c,
                           y_c=y_c,
                           z_c=z_c,
                           lvl_c=lvl_c,
                           vx_c=vx_c,
                           vy_c=vy_c,
                           vz_c=vz_c,
                           m_c=m_c,
                           T_c=T_c,
                           metal_c=metal_c,
        
                           pos_ctr=pos_ctr,
                           vel_ctr=vel_ctr,
        
                           l_simbox=param.l_simbox,
                           rotating_gal=param.rotating_gal,
                           R_upp=param.R_upp,
                           R_low=param.R_low,
                           A_upp=param.A_upp,
                           mass_weight=param.mass_weight,
                           par_mode=param.par_mode,
                           sig_s = param.sig_s, #[km/s]
                           old_fbase = param.repo,
                           old_fname = param.opart_name,
                           add_young=param.add_young,
                           young_fbase = param.repo,
                           young_fname = param.ypart_name,
                           Age_Cut = param.Age_Cut, #[Myr]
                           SED_type_y=param.SED_type_y,
                           import_velocity=param.import_velocity,
                           import_dispersion=param.import_dispersion,
                           gas_fbase = param.repo,
                           gas_fname = param.gas_name,
                           fact_dx = param.fact_dx,
                           sig_c = param.sig_c, #[km/s]
                           metal_dependent_dust = param.metal_dependent_dust,
                           Zsun = param.Zsun,
                           mtd_mod = param.mtd_mod,
                           import_velocity_gas=param.import_velocity_gas,
                           import_dispersion_gas=param.import_dispersion_gas,
                           cell_medium=param.cell_medium,
                           dx_s = param.dx_s, #[pc]
                           adaptive_smoothing=param.adaptive_smoothing,
                           n_targ_nbr=param.n_targ_nbr,
                           n_job_smoothing=param.n_job_smoothing
                          )
    
    make_SKI.writing_ski(
                   mode=param.mode,
                   N_phot=param.N_phot,
                   rnd_seed=param.rnd_seed,
                   wavelengthOutputStyle=param.wavelengthOutputStyle,
                   fluxOutputStyle=param.fluxOutputStyle,
                   z_red=param.z_red,
                   H_0=param.H_0,
                   Om_m=param.Om_m,
                   source_unit=param.source_unit,
                   wv_source_min=param.wv_source_min, 
                   wv_source_max=param.wv_source_max, 
                   wv_source_del=param.wv_source_del,
                   source_bias=param.source_bias, 
                   bias_type=param.bias_type,
                   bias_unit=param.bias_unit,
                   wv_bias_min=param.wv_bias_min,
                   wv_bias_max=param.wv_bias_max,
                   import_velocity=param.import_velocity,
                   import_dispersion=param.import_dispersion,
                   source_weight=param.source_weight,
                   bias_deg=param.bias_deg,
                   smoothing_kernel_o=param.smoothing_kernel_o,
                   SED_type_o=param.SED_type_o,
                   SED_IMF_o=param.SED_IMF_o,
                   SED_res_o=param.SED_res_o,
                   add_young=param.add_young,
                   smoothing_kernel_y=param.smoothing_kernel_y,
                   SED_type_y=param.SED_type_y,
                   SED_IMF_y=param.SED_IMF_y,
                   SED_res_y=param.SED_res_y,
                   forceScattering=param.forceScattering,
                   minWeightReduction=param.minWeightReduction,
                   minScattEvents=param.minScattEvents,
                   pathLengthBias=param.pathLengthBias,
                   storeRadiationField=param.storeRadiationField,
                   radfield_type=param.radfield_type,
                   radfield_unit=param.radfield_unit,
                   radfield_min_wv=param.radfield_min_wv, 
                   radfield_max_wv=param.radfield_max_wv, 
                   radfield_N_wv=param.radfield_N_wv, 
                   med_type=param.med_type,
                   med_Mfrac=param.med_Mfrac,
                   importMetallicity=param.importMetallicity,
                   importTemperature=param.importTemperature,
                   maxTemperature=param.maxTemperature, 
                   importVelocity_med=param.import_velocity_gas,
                   importMagneticField=param.importMagneticField,
                   importVariableMixParams=param.importVariableMixParams,
                   smoothing_kernel_med=param.smoothing_kernel_med, 
                   Dust_Type=param.Dust_Type,
                   Dust_Env=param.Dust_Env,
                   N_Si=param.N_Si,
                   N_C=param.N_C,
                   N_PAH=param.N_PAH,
                   numDensitySamples=param.numDensitySamples,
                   numPropertySamples=param.numPropertySamples,
                   aggregateVelocity=param.aggregateVelocity,
                   grid_unit=param.grid_unit,
                   grid_minX=param.grid_minX,
                   grid_maxX=param.grid_maxX,
                   grid_minY=param.grid_minY,
                   grid_maxY=param.grid_maxY,
                   grid_minZ=param.grid_minZ,
                   grid_maxZ=param.grid_maxZ,
                   Tree_Type=param.Tree_Type,
                   Tree_Policy=param.Tree_Policy,
                   minLevel=param.minLevel,
                   maxLevel=param.maxLevel,
                   numExtraLevels=param.numExtraLevels,
                   maxDustFraction=param.maxDustFraction,
                   maxDustOpticalDepth=param.maxDustOpticalDepth,
                   grid_wv_test_unit=param.grid_wv_test_unit,
                   grid_wv_test=param.grid_wv_test,
                   maxDustDensityDispersion=param.maxDustDensityDispersion,
                   maxElectronFraction=param.maxElectronFraction,
                   maxGasFraction=param.maxGasFraction,
                   def_wv_grid=param.def_wv_grid,
                   includeGALEX=param.includeGALEX,
                   includeSDSS=param.includeSDSS,
                   include2MASS=param.include2MASS,
                   includeWISE=param.includeWISE,
                   includeHERSCHEL=param.includeHERSCHEL,
                   inst_dist_unit=param.inst_dist_unit,
                   inst_dist=param.inst_dist,
                   inc_default=param.inc_default,
                   azm_default=param.azm_default,
                   roll_default=param.roll_default,
                   inc_min=param.inc_min,
                   inc_max=param.inc_max,
                   inc_del=param.inc_del,
                   azm_min=param.azm_min,
                   azm_max=param.azm_max,
                   azm_del=param.azm_del,
                   recordComponents=param.recordComponents,
                   recordPolarization=param.recordPolarization,
                   recordStatistics=param.recordStatistics,
                   numScatteringLevels=param.numScatteringLevels,
                   fov_unit=param.fov_unit,
                   fov_X=param.fov_X,
                   fov_Y=param.fov_Y,
                   pscale_X=param.pscale_X, 
                   pscale_Y=param.pscale_Y,
                   centre_X=param.centre_X,
                   centre_Y=param.centre_Y,
                   inst_2d_sed_on=param.inst_2d_sed_on,
                   inst_2d_sed_unit=param.inst_2d_sed_unit,
                   inst_2d_sed_type = param.inst_2d_sed_type,
                   inst_2d_sed_min_wv = param.inst_2d_sed_min_wv,
                   inst_2d_sed_max_wv = param.inst_2d_sed_max_wv,
                   inst_2d_sed_N_wv = param.inst_2d_sed_N_wv,
                   inst_2d_sed_min_wv_sub = param.inst_2d_sed_min_wv_sub,
                   inst_2d_sed_max_wv_sub = param.inst_2d_sed_max_wv_sub,
                   inst_2d_sed_N_wv_sub = param.inst_2d_sed_N_wv_sub,
                   inst_2d_sed_repo=param.inst_2d_sed_repo,
                   inst_2d_sed_fname=param.inst_2d_sed_fname,
                   inst_2d_sed_relHW=param.inst_2d_sed_relHW,
                   save_1d_sed=param.save_1d_sed,
                   aperture_unit=param.aperture_unit, 
                   aperture_min=param.aperture_min, 
                   aperture_max=param.aperture_max,
                   aper_del=param.aper_del,
                   inst_sed_grid_type=param.inst_sed_grid_type,
                   inst_sed_grid_unit=param.inst_sed_grid_unit,
                   inst_sed_min_wv=param.inst_sed_min_wv,
                   inst_sed_max_wv=param.inst_sed_max_wv,
                   inst_sed_N_wv=param.inst_sed_N_wv,
                   inst_sed_min_wv_sub=param.inst_sed_min_wv_sub,
                   inst_sed_max_wv_sub=param.inst_sed_max_wv_sub,
                   inst_sed_N_wv_sub=param.inst_sed_N_wv_sub,
                   repo=param.repo,
                   opart_name=param.opart_name,
                   ypart_name=param.ypart_name,
                   gas_name=param.gas_name,
                   ski_name=param.ski_name,
                   indent=param.indent
               )

    
    if param.exe_skirt == True:
        execute(repo=param.repo,
                repo_output=param.repo_output,
                ski_name=param.ski_name,
                skirt_dir=param.skirt_dir,
                N_thread=param.N_thread)
    
    return 

















