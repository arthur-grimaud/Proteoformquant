# """Store the deleted methods"""

# # ---------------------------------- Outlier detection --------------------------------- #

#     def exclude_outlier_non_significant(self, auc_percent_threshold):
#         "Exclude datapoint on the tails of the elution profile model that would represent less than X%  of the area under the curve"

#         psms_included = self.psms_included
#         psms_outliers = self.psms_outliers

#         for psm in self.psms_included:

#             auc_l_r = (
#                 self.get_auc(-np.inf, psm.spectrum.get_rt()),
#                 self.get_auc(psm.spectrum.get_rt(), np.inf),
#             )
#             auc_tot = sum(auc_l_r)
#             if auc_tot > 0:
#                 auc_min = min(auc_l_r)

#                 auc_percent = (auc_min / auc_tot) * 100

#                 if auc_percent < auc_percent_threshold / 2:  # diveded by 2 because two tailed
#                     psms_outliers.append(psm)

#         psms_included = [psm for psm in self.psms_included if psm not in psms_outliers]

#         return psms_outliers, psms_included

#     def exclude_outliers_mean_method(self):
#         "Try imrove the fit of the curve by iteratively removing datapoints the furthest from the RT mean"

#         subsets_psms = [self.psms_included]
#         subsets_rt = [[psm.spectrum.get_rt() for psm in self.psms_included]]
#         subsets_index = [0]

#         # print(subsets_psms[0][0].get_modification_brno())

#         # Create a list of psms subsets
#         for i in range(0, len(subsets_psms[0]) - 5):
#             m = mean(subsets_rt[i])
#             outIndex = subsets_rt[i].index(max(subsets_rt[i], key=lambda x: abs(x - m)))
#             subsets_psms.append([v for y, v in enumerate(subsets_psms[i]) if y != outIndex])
#             subsets_rt.append([v for y, v in enumerate(subsets_rt[i]) if y != outIndex])
#             subsets_index.append(i + 1)

#         # Compute fit score for each subsets
#         subsets_scores = []
#         for psmsSubset in subsets_psms:
#             data_yT = np.array([psm.spectrum.get_rt() for psm in psmsSubset])
#             yDataT = np.array([psm.get_prec_intens_ratio() for psm in psmsSubset])
#             # refit the curve to subset
#             param_estimated, param_fitted, score_estimated, score_fitted = self.fit_skew_normal(
#                 data_yT, yDataT
#             )
#             # store score and subset

#             subsets_scores.append(score_fitted)

#         # Find best score that retain the maximum number of psms
#         kn = KneeLocator(
#             subsets_index,
#             subsets_scores,
#             S=2,
#             curve="concave",
#             direction="increasing",
#             interp_method="polynomial",
#             polynomial_degree=2,
#         )
#         index = kn.knee
#         # print(index)
#         # Return best fit results
#         if index != None and index != 0:
#             try:
#                 data_x = np.array([psm.spectrum.get_rt() for psm in subsets_psms[index]])
#                 data_y = np.array([psm.get_prec_intens_ratio() for psm in subsets_psms[index]])
#                 param_estimated, param_fitted, score_estimated, score_fitted = self.fit_skew_normal(
#                     data_x, data_y
#                 )
#                 psms_outliers = [
#                     psm for psm in self.psms_included if psm not in subsets_psms[index]
#                 ]  # TODO check whether outlier needs to be removed for proteoform.linkedPSMSs
#                 psms_included = [psm for psm in self.psms_included if psm in subsets_psms[index]]
#             except (TypeError):
#                 print("TYPE ERROR: could not optimize fit by removing data points")
#                 psms_outliers = []
#         else:
#             psms_outliers = []
#             psms_included = self.psms_included
#             print("NO OPTI: could not optimize fit by removing data points")
#             pass

#         # print(subsets_index)
#         # print(subsets_scores)
#         # print(score_fitted)
#         return param_estimated, param_fitted, score_estimated, score_fitted, psms_outliers, psms_included


#     def __skewnormal_cdf(self, x, m, s, a, k, range_start, range_end):
#         values = []
#         for value in x:
#             integral = quad(lambda x: self.skewnormal(x, m, s, a, k), range_start, value)[0]
#             normalized = (
#                 integral / quad(lambda x: self.skewnormal(x, m, s, a, k), range_start, range_end)[0]
#             )
#             values.append(normalized)
#         return np.array(values)


#             # def __estimate_initial_parameters(self, model, data_x, data_y):

#     #     result = differential_evolution(self.__MSPD, self.__get_parameter_bounds_skewnormal(data_x, data_y), args =(model, data_x, data_y), seed=1)
#     #     return result.x


# f = open(f"quant_initial_{output_prefix}.csv", "w")
#     for proteoform in run.proteoforms.values():
#         if proteoform.update_proteoform_total_intens(method="precursor") > 0:
#             f.write(
#                 "".join(
#                     [
#                         str(proteoform.get_modification_brno()),
#                         ",",
#                         str(proteoform.get_peptide_sequence()),
#                         ",",
#                         str(proteoform.update_proteoform_total_intens(method="precursor")),
#                         ",",
#                         str(proteoform.update_proteoform_total_intens(method="AUC")),
#                         ",",
#                         str(proteoform.is_ambiguous()),
#                         "\n",
#                     ]
#                 )
#             )
#     f.close()
