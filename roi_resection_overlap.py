#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Created on Tue Jul 23 11:32:07 2024
@authors: Leela Srinivasan, Hera Akmal

Prerequisites: AFNI/SUMA, Freesurfer, generated resection mask
Description: Calculates overlap between resection and regions of interest. 
             Currently calculating overlap on the piriform cortex, hippocampus, and amygdala.
Usage: ./roi_resection_overlap [p number] [lh/rh (resection lateralization)]

"""

# Import processing modules
import sys
import os
import subprocess

# Read in command line arguments
subj=sys.argv[1]
rsxn_hemi=sys.argv[2]

# Define paths
working_dir = 'insert_working_directory'
fs_dir = 'freesurfer_directory'
subj_fs_dir = 'subject_freesurfer_directory'
suma_dir = 'subject_freesurfer_SUMA_directory'
session='subject_freesurfer_session'

# Locate SUMA surface volume / preoperative t1 / resection mask
surfvol_fname ='surface_volume_nifti_file'
t1 = 'preoperative_t1_nifti_file'
rsxn_msk ='resection_mask_nifti'

os.chdir(working_dir)
for hemi in ['lh', 'rh']:
    
    #Step 1: Map HCP-MMP1 atlas onto patient MRI (downloaded from atlas page and moved into fsaverage folder; 
    #use different atlas if desired)
    sval='{}/fsaverage/label/{}.HCP-MMP1.annot'.format(fs_dir, hemi)
    tval='{}/sub-{}_{}/label/{}.HCP-MMP1.annot'.format(fs_dir, subj, session, hemi)
    
    cmd1="mri_surf2surf --srcsubject fsaverage --trgsubject sub-{}_{} --hemi {} --sval-annot {} --tval {}"
    cmd1=cmd1.format(subj, session, hemi, sval, tval)
    subprocess.run(cmd1,shell=True)
    
    #Step 2: Convert annotation files to GIFTI files 
    white_surf = os.path.join(subj_fs_dir, 'surf', '{}.white'.format(hemi))
    gifti_output = os.path.join(working_dir, '{}.HCP-MMP1.gii'.format(hemi))

    cmd2 = "mris_convert --annot {} {} {}"
    cmd2=cmd2.format(tval, white_surf, gifti_output)
    subprocess.run(cmd2, shell=True)

#Step 3: Convert anatomical parcellation from the cortical surface (annotation file) to volumetric segmentation
#hcp_mmp1_vol = os.path.join(subj_fs_dir, 'mri', 'HCP-MMP1_vol.mgz')
hcp_mmp1_vol = 'HCP-MMP1_vol.mgz'
cmd3 = "mri_aparc2aseg --s sub-{}_{} --annot HCP-MMP1 --o {}"
cmd3=cmd3.format(subj, session, hcp_mmp1_vol)
subprocess.run(cmd3, shell=True)

#Step 4: Convert .mgz file to .nii
hcp_mmp1_vol_nii = 'HCP-MMP1_vol.nii'
cmd4 = "mri_convert {} {}"
cmd4=cmd4.format(hcp_mmp1_vol, hcp_mmp1_vol_nii)
subprocess.run(cmd4, shell=True)

#Step 5: Extract volumes for piriform cortex, hippocampus, and amygdala
pir_stats = 'piriform_cortex.stats'
hip_amy_stats = 'hippocampus_amygdala.stats'

cmd5a = "mri_segstats --seg {} --id 1110 2110 --sum {}"
cmd5a=cmd5a.format(hcp_mmp1_vol, pir_stats)
subprocess.run(cmd5a, shell=True)

cmd5b = "mri_segstats --seg {} --id 17 53 18 54 --sum {}"
cmd5b=cmd5b.format(hcp_mmp1_vol, hip_amy_stats)
subprocess.run(cmd5b, shell=True)

#Step 6: Create binary masks for piriform cortex, hippocampus, and amygdala
#Mask for piriform is created using HCP-MMP1 volume whereas hippocampus and amygdala are created using default 
#aseg volume; change as needed
pir_mask = '{}_piriform_mask.mgz'.format(rsxn_hemi)
hip_mask = '{}_hippocampus_mask.mgz'.format(rsxn_hemi)
amy_mask = '{}_amygdala_mask.mgz'.format(rsxn_hemi)
aseg_vol = os.path.join(subj_fs_dir, 'mri', 'aseg.mgz')

cmd6a = "mri_binarize --i {} --match {} --o {}"
cmd6a = cmd6a.format(hcp_mmp1_vol, '1110' if rsxn_hemi == 'lh' else '2110', pir_mask)
subprocess.run(cmd6a, shell=True)

cmd6b = "mri_binarize --i {} --match {} --o {}"
cmd6b = cmd6b.format(aseg_vol, '17' if rsxn_hemi == 'lh' else '53', hip_mask)
subprocess.run(cmd6b, shell=True)

cmd6c = "mri_binarize --i {} --match {} --o {}"
cmd6c = cmd6c.format(aseg_vol, '18' if rsxn_hemi == 'lh' else '54', amy_mask)
subprocess.run(cmd6c, shell=True)

#Step 7: Extract volume for resection mask
rsxn_msk_mgz = 'rsxn.msk.mgz'
rsxn_msk_stats = 'rsxn_msk_stats.txt'

cmd7a = "mri_convert {} {}"
cmd7a=cmd7a.format(rsxn_msk, rsxn_msk_mgz)
subprocess.run(cmd7a, shell=True)

cmd7b = "mri_segstats --seg {} --sum {}"
cmd7b=cmd7b.format(rsxn_msk_mgz, rsxn_msk_stats)
subprocess.run(cmd7b, shell=True)

#Step 8: Convert .mgz files to .nii format
pir_mask_nii = '{}_piriform_mask.nii'.format(rsxn_hemi)
hip_mask_nii = '{}_hippocampus_mask.nii'.format(rsxn_hemi)
amy_mask_nii = '{}_amygdala_mask.nii'.format(rsxn_hemi)

cmd8a = "mri_convert {} {}"
cmd8a=cmd8a.format(pir_mask, pir_mask_nii)
subprocess.run(cmd8a, shell=True)

cmd8b = "mri_convert {} {}"
cmd8b=cmd8b.format(hip_mask, hip_mask_nii)
subprocess.run(cmd8b, shell=True)

cmd8c = "mri_convert {} {}"
cmd8c=cmd8c.format(amy_mask, amy_mask_nii)
subprocess.run(cmd8c, shell=True)

#Step 9: Align resection mask
cmd9a="3dAllineate -base {} -source sub-{}_{}_SurfVol.nii -prefix aligned+orig -1Dmatrix_save fs_to_anat"
cmd9a=cmd9a.format(t1, subj,session)
subprocess.run(cmd9a,shell=True)

#Iterate through ROI mask nifti files to apply transformation and generate aligned masks
for mask_nii in [pir_mask_nii, hip_mask_nii, amy_mask_nii]:

    # Apply transformation matrix to the ROI mask
    cmd9b="3dAllineate -base {} -source {} -1Dmatrix_apply fs_to_anat.aff12.1D -prefix tmp+orig"
    cmd9b=cmd9b.format(t1,mask_nii)
    subprocess.run(cmd9b,shell=True)
        
    # Verify allineation
    cmd9c="3dcalc -a tmp+orig -expr 'ispositive(a-0.1)' -prefix tmp2.nii"
    subprocess.run(cmd9c,shell=True)
            
    # Fill holes
    cmd9d="3dmask_tool -input tmp2.nii -prefix al_{} -fill_holes"
    cmd9d=cmd9d.format(mask_nii)
    subprocess.run(cmd9d,shell=True)
    
    # Delete temporary files
    for file in os.listdir(os.getcwd()):
        if 'tmp+orig' in file or 'tmp2' in file:
            os.remove(file)

#Step 10: Compute overlap and output to text file
overlap_output = os.path.join(working_dir, 'overlap_output.txt')

#Iterate through ROIs to calculate overlap of each with resection mask
with open(overlap_output, 'w') as f:
    for roi in ['piriform', 'hippocampus', 'amygdala']:
        cmd10 = "3dOverlap -save {}_{}_overlap al_{}_{}_mask.nii rsxn.msk.nii"
        cmd10=cmd10.format(rsxn_hemi, roi, rsxn_hemi, roi)
        result10 = subprocess.run(cmd10, shell=True, capture_output=True, text=True)
        
        # Write label and the result to the text file
        f.write("Number of voxels overlapping with {} {}:\n".format(rsxn_hemi, roi))
        f.write(result10.stdout)
        f.write('\n')
        
#Step 11: Compute number of voxels in masks and output to text file
num_voxels_output = os.path.join(working_dir, 'num_voxels_output.txt')

#Extract number of voxels in resection mask
with open(num_voxels_output, 'w') as f:
    cmd11a = "3dBrickStat -count -non-zero {}"
    cmd11a=cmd11a.format(rsxn_msk)
    result11a = subprocess.run(cmd11a, shell=True, capture_output=True, text=True)
    
    # Write label and the result to the text file
    f.write("Number of voxels in resection mask:\n")
    f.write(result11a.stdout)
    f.write('\n')
    
    #Iterate through ROIs to extract number of voxels in each
    for roi in ['piriform', 'hippocampus', 'amygdala']:
        cmd11b = "3dBrickStat -count -non-zero al_{}_{}_mask.nii"
        cmd11b=cmd11b.format(rsxn_hemi, roi)
        result11b = subprocess.run(cmd11b, shell=True, capture_output=True, text=True)
        
        # Write label and the result to the text file
        f.write("Number of voxels in {} {}:\n".format(rsxn_hemi, roi))
        f.write(result11b.stdout)
        f.write('\n')
        