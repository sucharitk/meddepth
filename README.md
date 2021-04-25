# meddepth

processing pipeline and statistical analysis code for the project on subjective experience of meditation depth (Meditation depth questionnaire) in long-term meditators and matched control participants along with EEG, heart and respiration data

preprint available at:
Katyal, S., & Goldin, P. (2020, October 12). Alpha and theta brain oscillations play a complementary role in inducing levels of altered consciousness during meditation practice. PsyArXiv. https://doi.org/10.31234/osf.io/hezfd


directory structure

'processing' contains matlab code for processing the minimally pre-processed data. the main script is meddepth_analysis.m, which calls all the other files in the directory 

'stats' contains the R code for doing the analysis and figures in the paper 

'helper' contains additional .m files called by the matlab code

supplementary data including audios for guided meditations and the baseline block used in the study are available at https://osf.io/sfkte/
