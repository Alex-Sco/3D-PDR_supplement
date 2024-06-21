# 3D-PDR_supplement

## read_chemistry.py

### Requirements:

1. Python3

2. The following Python3 packages:
* numpy 
* pandas
* matplotlib 
* os 
* sys 

3. A directory named 'chemical_colors', aligned with the colorpwd in the procedure (this directory name can be modified as needed).

4. [modelname].chemistry.fin files for the models you wish to analyze.

### How to use & Examples:

There are two ways to use it: 
I. Input parameters externally:

In python3 environment: 
```
python3 read_chemistry.py [Species name] [model1 name] [model2 name] [model3 name] [Figure name] 
```

or, in the Ipython3 environment:
```
In [1]: run read_chemistry.py [Species name] [model1 name] [model2 name] [model3 name] [Figure name]
```

For example, if we want to compare three model results named "O0C0D0_15.chemistry.fin" "O0C1D0_15.chemistry.fin" and "O0CupD0_15.chemistry.fin" in the local directory, we can:

Step 1. Create a directory named 'chemical_colors' (here in the local directory)

```
mkdir chemical_colors
```

Step 2. Check the model path and the color file path in the script (go into the script and modify the path if you need)

``` 
colorpwd="./chemical_colors/"   #the directory where the color-reactions file is stored.
modelpwd="./" #the directory where the model results are stored. 
```

Step 3. Run the script (here plot reactions of CO as an example)

```
python3 read_chemistry.py "CO" "O0C0D0_15" "O0C1D0_15" "O0CupD0_15" "O0models"
```

Then we will get:
* Two files named 'For_reaction_colors_CO.dat' and 'Des_reaction_colors_CO.dat' in the chemical_colors directory, store the information on CO formation/destruction reactions and their corresponding colors.
* One file named 'CO_O0models.pdf' in the local directory, shows the comparison of CO formation and destruction reactions between the three models.


### Current limitations:

1. To ensure visibility against a white background and distinction from other colors, the current version includes only 21 colors. This means that for each species, the number of individual formation or destruction reactions displayed is up to 21.

2. While the code generates a file storing the colors and their corresponding reactions for each species, these files are species-specific. This ensures that within a species, a reaction is consistently represented by one color. However, the same reaction may be represented by different colors across figures of different species. This is what I planned to improve.

3. The current version must compare three models at a time. This can be modified by changing the relevant plotting section in the code.


******
If you have any questions, please contact me: yichen_sun@smail.nju.edu.cn

Many thanks for your support! :)
