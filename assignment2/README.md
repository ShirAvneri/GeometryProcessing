# Assignment 2
  
## Required results

### Mandatory Tasks
1) Show the visualization of the constrained points for the 'cat.off' point cloud.  
![Cat Constrained Points](images/CatConstrained.png)  
![Cat Constrained Points](images/CatConstrained2.png)  

2) Show screenshots of the grid with nodes colored according to their implicit function values (cat.off and luigi.off).  
#### cat.off
![Cat Grid](images/CatGrid.png)  
#### luigi.off
![Luigi Grid](images/LuigiGrid.png)  

3) Show screenshots of the reconstructed surfaces. Experiment with different parameter settings: grid resolution (also anisotropic in the 3 axes), Wendland function radius, polynomial degree. Add all these settings to the GUI to ease experimentation. Briefly summarize your observations and save the reconstructed models in the off format for every point-cloud dataset provided (assignment2/results).
### Key Observations  
* Resolution - the higher the resolution, the better the reconstructed models. However, if it is not uniform, the models twitch.
* Wendland Radius - a high Wendland Radius can lead to parts of the model disappearing, and a low one can lead to holes and distortion in the model.
* Polynomial Degrees - 1 and 2 are likely to produce artifacts around the model.

### Experimenting with the Cat Model  
#### Resolution: 10x10x10, Wendland Radius: 0.100, Polynomial Degree: 0  
![Reconstructed Cat 1](images/ReconstructedCat1.png)  
#### Resolution: 20x20x20, Wendland Radius: 0.100, Polynomial Degree: 0  
![Reconstructed Cat 2](images/ReconstructedCat2.png)  
#### Resolution: 30x30x30, Wendland Radius: 0.100, Polynomial Degree: 0  
![Reconstructed Cat 3](images/ReconstructedCat3.png)  
#### Resolution: 20x20x20, Wendland Radius: 0.200, Polynomial Degree: 0  
![Reconstructed Cat 4](images/ReconstructedCat4.png)  
#### Resolution: 20x20x20, Wendland Radius: 0.050, Polynomial Degree: 0  
![Reconstructed Cat 5](images/ReconstructedCat5.png)  
#### Resolution: 20x20x20, Wendland Radius: 0.100, Polynomial Degree: 1  
![Reconstructed Cat 6](images/ReconstructedCat6.png)  
#### Resolution: 20x20x20, Wendland Radius: 0.100, Polynomial Degree: 2  
![Reconstructed Cat 7](images/ReconstructedCat7.png)  
#### Resolution: 25x20x20, Wendland Radius: 0.100, Polynomial Degree: 0  
![Reconstructed Cat 8](images/ReconstructedCat8.png)  

### Experimenting with Other Models  
#### Resolution: 20x20x20, Wendland Radius: 0.100, Polynomial Degree: 0  
![Reconstructed Bunny-1000](images/ReconstructedBunny1000.png)  
#### Resolution: 37x37x37, Wendland Radius: 0.100, Polynomial Degree: 0  
![Reconstructed Luigi](images/ReconstructedLuigi.png)  
#### Resolution: 20x20x20, Wendland Radius: 0.200, Polynomial Degree: 1  
![Reconstructed Sphere](images/ReconstructedSphere.png)  

4) Theory question: Save your notes to assignment2/results and add a link to this page.  
[Solution](results/Theory_Question.pdf)
