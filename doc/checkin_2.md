# Check-In 2
## Addressing Prior Feedback
- How will you test and document tests in this repo about the rewritting of PerturbSeq?
- The original PerturbSeq documentation contains ipynb notebooks with sample usage of the package. Additionally, there is a separate repo that is self-contained and has all the data required for a perturbseq demo, which we think is the best initial test: https://github.com/thomasmaxwellnorman/perturbseq_demo (i.e., we will write the perturbseq_demo Python notebook into an R notebook and test to see if the demo can be recreated)

- ***add comment about cell ranger data***
## Progress
Dan Le:
- 

Jander Kugelman:
- Wrote transformers.R without tsne implementation
- About halfway through cell_population.R rewrite. This will likely be the most problematic package to rewrite and test because of how long it is.
## Project Organization
- figures we are attempting to recreate have been moved to /images
- progress on perturbseq package is in /src/perturbseqR/R
## Struggles/Questions
Jander: 
- I am still in the process of rewriting the last script (src/perturbseqR/R/cell_population.R). The only struggle I have right now with it is that it is about 1000 lines of code, so it is taking me awhile to make sure I am rewriting all of the functions correctly.
- Dan figured out that the tsne package for python used in the transformers.py script in the perturbseq package is no longer functional, so I only wrote the sklearn