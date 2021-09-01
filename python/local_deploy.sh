#!/bin/bash 

location=$1
cp plotCaptureHeatmapQC.v-11.py plotCaptureHeatmapQC.cc.py
cp combineCorsivWIRITC.v-1.py combineCorsivWIRITC.cc.py
cp recodeCoRSIV.v-1.py recodeCoRSIV.cc.py
cp captureBetaVsGeneExpression.v-7.py captureBetaVsGeneExpression.cc.py
cp corsivCaptureVersionsCheck.v-1.py corsivCaptureVersionsCheck.cc.py
cp plotCaptureDepth.v-4.py plotCaptureDepth.cc.py
cp summarizeGeneBetaCorr.v-6.py summarizeGeneBetaCorr.cc.py
cp submiCaptureQC.v-1.py submitCaptureQC.cc.py
cp submitCaptureQCRatio.v-1.py submitCaptureQCRatio.cc.py
cp captureQC.v-1.py captureQC.cc.py 
cp captureQCRatio.v-1.py captureQCRatio.cc.py 
cp featureDensityPlot.v-4.py featureDensityPlot.cc.py
cp driver_eMatrixQTL.v-1.py driver_eMatrixQTL.cc.R
cp splitCorsivsByChrom.v-1.py splitCorsivsByChrom.cc.py
cp runCoRSIVmQTL.v-3.py runCoRSIVmQTL.cc.py
cp submitEMatrixQTL.v-2.py submitEMatrixQTL.cc.py
cp annotate_mQTL_distance.v-3.py annotate_mQTL_distance.cc.py
cp annotateSNVCorsivFDR.v-3.py annotateSNVCorsivFDR.cc.py
cp simes_correction.v-1.R simes_correction.cc.R
cp addBetaToSimes.v-1.py addBetaToSimes.cc.py
cp addDistanceToSimes.v-1.py addDistanceToSimes.cc.py
cp select_SNP_CoRSIV_by_FDR.v-1.py select_SNP_CoRSIV_by_FDR.cc.py
cp distance_histogram_CoRSIV_SNP.v-1.py distance_histogram_CoRSIV_SNP.cc.py
cp annotate_corsiv_snp_distance.v-1.py annotate_corsiv_snp_distance.cc.py
cp hapblockCoRSIV.v-7.py hapblockCoRSIV.cc.py
cp submitHapblockCoRSIV.v-1.py submitHapblockCoRSIV.cc.py
cp glm_corsiv_hapsum.v-2.R glm_corsiv_hapsum.cc.R
cp glm_corsiv_allele.v-3.R glm_corsiv_allele.cc.R
cp subsetRepeats.v-1.py subsetRepeats.cc.py
cp repeatDensityVsCtrl.v-4.py repeatDensityVsCtrl.cc.py
cp subset_corsiv_repeat_odds_ratio.v-1.py subset_corsiv_repeat_odds_ratio.cc.py
cp annotateCorsivIndex.v-1.py annotateCorsivIndex.cc.py
cp combineBestMQTLCoRSIV.v-3.py combineBestMQTLCoRSIV.cc.py
cp summarizeMQTLStats.v-2.py summarizeMQTLStats.cc.py
cp summarizeSelectionScore.v-1.py summarizeSelectionScore.cc.py
cp filterGTExSNPs.v-3.py filterGTExSNPs.cc.py
cp manhattan_mQTL_plot.traditional.v-4.R manhattan_mQTL_plot.traditional.cc.R
cp select_mQTL_Quartiles.v-1.py select_mQTL_Quartiles.cc.py
cp nearestCoRSIVGene.v-1.py nearestCoRSIVGene.cc.py
cp medianMQTLCandidateDistance.v-1.py medianMQTLCandidateDistance.cc.py

chmod +x *.cc.py *.cc.R
cp -f *.cc.py *.cc.R ${location}/
# rm -f *.cc.py *.cc.R
