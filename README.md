# Data & analysis

## Vogels et al., 2019. 'Unsolved Mystery: Arbovirus co-infection and co-transmission: a neglected public health concern?'. PLOS Biology

### Overlapping Zika, dengue, and chikungunya virus outbreaks in the Americas (Figure 1)

Data was obtained from the Pan-American Health Organization. Dengue and Chikungunya cases were extracted from the cummulative case reports, and thus may not accurately reflect the true number of reported cases per month. Zika cases were extracted from bar graphs showing cases per week. See [this blog post](http://andersen-lab.com/paho-zika-cases/) for info on how monthly Zika cases were extracted.

Files:
* CHIKV_2016-2017.xlsx
* DENV_2016-2017.xlsx
* ZIKV_2016-2017.xlsx


### Summary of arbovirus human co-infection studies (Table 1)

We compiled data from reports of human Zika, dengue, and/or chikungunya virus co-infections over the last 50 years and the clinical presentation associated with them.

File:
* Human_co-infections.xlsx

### Impact of mosquito co-infection on transmission (Figure 3B)

We compiled data on mosquito transmission from studies that made a direct comparison between mosquitoes exposed to a single or multiple viruses. Transmission of co-exposed mosquitoes was calculated relative to single-exposed mosquitoes, with relative transmission being defined as: transmission rate of virus ‘X’ in mosquitoes co-exposed to virus X and Y / transmission rate virus X in single-exposed mosquitoes. Transmission is expressed as the percentage of mosquitoes with virus in their saliva out of the total number of exposed mosquitoes. Relative transmission of 1 indicates that no difference was observed between transmission rates of single-exposed or co-exposed mosquitoes.

File:
* Mosquito_co-transmission.xlsx

### Model-predicted prevalence over time of two sequentially invading arboviruses (Box)
We used a deterministic SIR-SI model to explore possible impacts that co-transmission from mosquitoes to humans may have on the overall dynamics of simultaneous arbovirus outbreaks. This model incorporates two viruses that have identical transmission parameters and recovery rates (for humans). The transmission parameters of the viruses are identical, and virus Y invades one month after virus X in a population of 1,000,000. Co-transmission from human to mosquito is fixed such that 60% of infectious bites on a co-infected human lead to co-transmission. Co-transmission from mosquito to human is varied between 0 and 100%. For further details see 'S2 File in the associated article (Vogels et al., 2019).

The model does not require any data inputs to run, but requires one input (Mosquito_To_Human.png) to generate the first figure, which should be in the current working folder. Model parameters can be varied near the top of the file. 	

Files:
* Co-infection_model_code.py
* Mosquito_To_Human.png


---
**Grubaugh Lab**  
Yale School of Public Health 
New Haven, CT, USA  
[grubaughlab@gmail.com](mailto:grubaughlab@gmail.com)
