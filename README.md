# FST_APP
FST_APP for Bioinformatics Course

This app is built to guide senior high school students through the process of determining the genetic distance among different populations of fish in Lake Mistassini in Eeyou Istchee. Hopefully, students will be able to use something like this to see the genetic data associated with the fishery that supports a large part of their community's economy. It would not stand alone, but be a part of a larger learning module.
Notes:
Some values are reactive even though this seems unnecessary at the moment (for example, the population (sites) checkboxes at the top of the plots tab). This is because at some point, other fish data sets will hopefully be added, and the sites may change, possibly to other nearby geographic locations. 

A minor bug in this App is that the number of fish designated as "kept" during the filtering process actually is increased when the filters are finalized. This is because the threshold that is set applied in the previous step is applied first during the finalization, thus restricting the number of loci available for each genotyping of each fish to fail at. Getting around this would require adding another aspect of reactivity that would likely cause more delay in computation times. And no doubt, that flopping fish is seriously annoying.

In addition to the sources that have been acknowledged in the code for this app (non-package functions and .gif found in server.r and ui.r files) much of the population genetics analysis pipeline was based on this helpful tutorial by Tom Jenkins: https://tomjenkins.netlify.app/2020/09/21/r-popgen-getting-started/.