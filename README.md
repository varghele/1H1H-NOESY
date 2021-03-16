# 1H1H-NOESY
## Script(s) to calculate 1H-1H-NOESY MAS NMR cross relaxation rates and plot them

If you use the python version, run the requirements.txt with ```pip install -r requirements.txt```.

An .exe version can be found on my GoogleDrive [here](https://drive.google.com/file/d/18Cr7hPeG7Gx46BN2Uoz7a0ObPgBQmqQO/view?usp=sharing) and was compiled for 32bit.

An examples can be found for the small molecule kinase inhibitor Ceritinib in the examples folder.

Script is fairly straightforward to use:
1. Integrate peaks in Topspin (you probably have done that already), and copy their values into the .txt files. For example into ```100ms.txt``` for 100ms mixing time.

2. Write the names of your peaks of the lipid and the ligand into ```ini_data\chemical.txt``` and ```ini_data\lipid.txt``` respectively, as well as the name of the lipid and ligand.

3. In ```plot_data\lipid_sorted.txt``` you can sort the peaks in the order you want. In ```plot_data\plot_data.txt``` you can set the bar colors.

4. Run the script.

5. In ```zresults\yourdrugname\...``` you can find the results. If you want to eliminate certain bars from the plot, set 1->0 in ```barplot_elimination.txt``` for the corresponding peak.

6. If you want to check the fits (done with Levenberg-Marquardt), and maybe take out a point there, simply navigate to ```x.xppm``` and modify the corresponding peak in  ```x.x_eliminator.txt``` from 1->0

7. Run script again, changes are applied.

8. If you want to do the plot yourself, the cross relaxation rates are in ```sigma.txt```, the errors in ```sigma_errors.txt```, and the correlation time in ```T_is.txt```

9. Good Luck.
-V
