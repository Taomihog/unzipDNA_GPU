## DNA unzipping curve prediction (GPU version )

### Instructions to Use (with Example Data)
1. The example genome is NEB H5alpha, the genbank file is in "examples" folder.
2. Run **parse_h5alpha.py** to parse the *200304_NEB_H5alpha_genbank.txt* and *200304_NEB_H5alpha.txt* to individual files. These files should be saved in folder "parsed"
3. Build the project on Windows. This .bat file will automatically run a test using the sequences in folder "parse"
```
> src/build.bat
```
4. To batch process sequence files, use this command
```
> main path/to/folder
```

### DNA unzipping curve prediction by GPU
(This is a learning project and my first time use CUDA 🧬)


### DNA unzipping theory

![pred_vs_exp](https://github.com/Taomihog/unzipDNA_GPU/blob/master/reference/DNA_unzipping_prediction_vs_sm_experiment.png)

**DNA unzipping experiment on a 4.4 kb DNA**. The program's prediction (black) aligns well with a single-molecule unzipping curve (blue).  
  
Further reading on DNA unzipping experiment and theory:  

[1] Essevaz-Roulet, Baptiste, Ulrich Bockelmann, and Francois Heslot. (1997) PNAS  
[2] Bockelmann, Ulrich, et al. (2002) Biophysical journal  
[3] Huguet, Josep M., et al. (2010) PNAS  
