# DNA Unzipping Curve Calculator (GPU version)  

This is my first time using CUDA. It runs pretty fast and can finish the whole *E.coli*  genome in about a minute on my GeForce RTX 3070 Laptop.  

## Instructions to Use (with Example Data)  

1. The example genome is NEB H5alpha, the genbank file is in "examples" folder.  
2. Run **parse_h5alpha.py** to parse the *200304_NEB_H5alpha_genbank.txt* and *200304_NEB_H5alpha.txt* to individual files. These files should be saved in folder "parsed"  
3. Build the project on Windows. This .bat file will automatically run a test using the sequences in "parse" folder  

```bash
> build.bat
```
To batch-process sequence files in your own folder, use this command  

```bash
> main path/to/your/own/folder
```

Result (I only showed the first 23 genes):  

![image](reference/result.svg)

Prediction vs Experiment:

![image](reference/theory_vs_experiment.png)

**DNA unzipping experiment vs theory**. The prediction (**${\color{black}-}$**) aligns well with the single-molecule unzipping curve (**${\color{red}-}$**).  

## Why CUDA  


Single-molecule approaches such as optical tweezers, can unzip a single dsDNA at single-molecular level.

![image](reference/sm_DNA_unzipping_exp_schematics.png)  

The theoretical prediction of an unzipping trace (*unzipping force vs total extension*) can be obtained from partition functin of the system. The partition function $Z$ at a total extension of the system $z$ is

$$Z(z) = \sum_{j=0}^{j_{max}}e^{-G(j,z)/kT}$$

where $G(j,z)$ is the total energy of the system. Force $F(z)$ can be obtained from partition function as follows:  

$$F(z) = -kT\frac{\partial }{\partial z}\mathrm{ln}Z$$

To calculate the force-extension curve, we need to obtain the energy $G(j,z)$ at every j and z. However, $G$ has a **non-analytical complex form**, while the scales of $j$ and $z$ are large (usually in a range of 1,000-10,000). Moreover, we need to calculate the unzipping traces for every gene in the genome. Even *E. Coli* has thousands of genes. Therefore, Calculation of $G$ for all $j$ and $z$ for all DNA sequences is better to be on a GPU.  

*(Besides using GPU, I also further sped up the calculation using a trick I developed in my previous project [unzipDNA_CPU](https://github.com/Taomihog/unzipDNA_CPU).)*

## Further Reading  

[1] Essevaz-Roulet, Baptiste, Ulrich Bockelmann, and Francois Heslot. (1997) PNAS  
[2] Bockelmann, Ulrich, et al. (2002) Biophysical journal  
[3] Huguet, Josep M., et al. (2010) PNAS  
