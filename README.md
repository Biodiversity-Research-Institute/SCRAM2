## What is SCRAM?

To estimate risk of avian collisions with offshore wind energy turbines in the U.S. Atlantic, a stochastic collision risk model using movement data from telemetry studies was developed for three federally protected bird species: the Roseate Tern (Sterna dougallii), Piping Plover (Charadrius melodus), and Red Knot (Calidris canutus). An online web application of the model, called Stochastic Collision Risk Assessment for Movement (SCRAM), and accompanying user manual have been made publicly available to help transparently estimate collision risk from offshore wind farms in the U.S. Atlantic.

 
## Study Background

Collision risk models are often used to estimate risk of avian collisions with offshore wind energy turbines. Such models typically use avian density data derived from observational survey datasets for specific locations, along with a suite of behavioral and site-specific variables that predict collision risk. However, very limited survey data are available for the Roseate Tern, Piping Plover, and Red Knot, all of which are federally protected under the U.S. Endangered Species Act. Thus, automated radio telemetry data and satellite telemetry data have been used to develop movement models that are linked to monthly population estimates, density estimates at specific flight heights, and other species- and site-specific characteristics (like number of turbines in a specified turbine array) to estimate collision risk for locations across a portion of the NES.

 
### The collision risk estimation process includes four main components:

(1) movement modeling to determine monthly occupancy rates over the NES,

(2) linking of monthly population estimates to occupancy rates to estimate species density across the NES,

(3) flight height estimation from Motus and GPS telemetry data to further refine the proportion of the population at risk of collision, and

(4) a collision risk model that uses density estimates at specific flight heights (along with a suite of other species- and location-specific parameters) to estimate collision for a specified turbine array.

SCRAM packages these model components into a web application for transparency and ease of use:

![SCRAM components](https://static.wixstatic.com/media/22a6e7_2274f8ca237f45baa0ea57c78d76607a~mv2.jpg/v1/crop/x_201,y_207,w_945,h_305/fill/w_944,h_305,al_c,q_80,enc_avif,quality_auto/SCRAM%20graphic%20for%20final%20report%20-%20draft%202.jpg)
 

## SCRAM Resources and Updates​

SCRAM 2

v2.1.8: Minor adjustments to Roseate Tern and Piping Plover Motus-based movement models after reassignment of 11 duplicate tag IDs for Roseate Terns and 10 duplicate tag IDs for Piping Plovers, and removal of 2 static Piping Plovers with long, uninformed data gaps

v2.1.6: Major update to improve movement models and add satellite tracking data to inform Red Knot movement/ flight height models


## Next Steps​

SCRAM updates will continue through at least 2026 with model improvements, bug fixes, and additional functionality. Updates for Phase 3 of SCRAM in 2025-2026 are expected to include:

Updated regional bird population size estimates provided by the U.S. Fish and Wildlife Service, based on the best available data

Improved precision of movement modeling estimates and updates to movement model structure

Incorporation of additional data into movement models, including new species, new tracking technologies, and new datasets for species and tracking technologies already included in SCRAM (as available)

Further development of the flight height modeling process as new data become available and methods are refined

 

## Acknowledgements

SCRAM was funded by the Bureau of Ocean Energy Management (BOEM) and the U.S. Fish and Wildlife Service (USFWS). Authors thank Wendy Walsh and Anne Hecht of the USFWS for their input on regional population size estimates for inclusion in the models, and David Bigger (BOEM) for his review of the models, user interface, and user manual. We would also like to thank five reviewers from the U.S. Geological Survey (USGS) for their input on an earlier version of SCRAM. We thank the funders and Principal Investigators of the red knot studies that collected satellite telemetry data used in SCRAM, including Stephanie Feign and Larry Niles (Wildlife Restoration Partnerships), Paul Smith (Environment and Climate Change Canada), Scott Lawton and Matt Overton (Dominion Energy, Coastal Virginia Offshore Wind Project), Candice Cook-Ohryn (Shell, Atlantic Shores Offshore Wind), and Kim Peters and Liz Gowell (Orsted, Ocean Wind). SCRAM builds on significant intellectual and analytical contributions from earlier collision risk models, including models by Band (2012), Masden (2015), and McGregor et al. (2018), as well as the current implementation of stochCRM in the stochLAB R package (Caneco et al. 2022).

​

## Disclaimers

Any use of trade, firm, or product names is for descriptive purposes only and does not imply endorsement by the U.S. Government. The findings and conclusions in SCRAM materials are those of the author(s) and do not necessarily represent the views of the U.S. Fish and Wildlife Service.
