## DNA Unzipping Curve Calculator (GPU version)  

This is a learning project and my first time using CUDA. It runs pretty fast and can finish the whole *E.coli*  genome in less than a minute.  

### Instructions to Use (with Example Data)  

1. The example genome is NEB H5alpha, the genbank file is in "examples" folder.  
2. Run **parse_h5alpha.py** to parse the *200304_NEB_H5alpha_genbank.txt* and *200304_NEB_H5alpha.txt* to individual files. These files should be saved in folder "parsed"  
3. Build the project on Windows. This .bat file will automatically run a test using the sequences in folder "parse"  

```bash
> src/build.bat
```

4. To batch-process sequence files in your own folder, use this command  

```bash
> main path/to/your/own/folder
```

### DNA unzipping curve prediction by GPU
(This is a learning project and my first time use CUDA 🧬)


### DNA unzipping theory

![pred_vs_exp](https://github.com/Taomihog/unzipDNA_GPU/blob/master/reference/DNA_unzipping_prediction_vs_sm_experiment.png)

**DNA unzipping experiment vs theory**. The prediction (black) aligns well with a single-molecule experimental data (blue).  


### DNA Unzipping Theory  


Further reading on DNA unzipping experiment and theory:  

[1] Essevaz-Roulet, Baptiste, Ulrich Bockelmann, and Francois Heslot. (1997) PNAS  
[2] Bockelmann, Ulrich, et al. (2002) Biophysical journal  
[3] Huguet, Josep M., et al. (2010) PNAS  
