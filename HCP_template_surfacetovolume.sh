export SUBJECTS_DIR="/.../recon"

mri_surf2surf --hemi lh \
              --srcsubject fsaverage \
              --sval-annot /.../lh.HCP2016infsaverage.annot \
              --trgsubject ch2 \
              --trgsurfval $SUBJECTS_DIR/ch2/label/lh.HCP.annot

mri_surf2surf --hemi rh \
              --srcsubject fsaverage \
              --sval-annot /.../rh.HCP2016infsaverage.annot \
              --trgsubject ch2 \
              --trgsurfval $SUBJECTS_DIR/ch2/label/rh.HCP.annot
              
tkregister2 --mov /.../recon/ch2/mri/orig.mgz \
            --noedit \
            --s ch2 \
            --regheader \
            --reg register_ch2.dat

mri_label2vol --annot /.../recon/ch2/label/lh.HCP.annot \
              --temp /.../recon/ch2/mri/orig.mgz \
              --subject ch2 \
              --hemi lh \
              --proj frac 0 1 0.1 \
              --fillthresh 0.3 \
              --o lh.HCP.ch2.nii \
              --reg register_ch2.dat

mri_label2vol --annot /.../recon/ch2/label/rh.HCP.annot \
              --temp /.../recon/ch2/mri/orig.mgz \
              --subject ch2 \
              --hemi rh \
              --proj frac 0 1 0.1 \
              --fillthresh 0.3 \
              --o rh.HCP.ch2.nii \
              --reg register_ch2.dat
