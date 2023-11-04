#rm(list = ls())
#mount bwFor directory to cluster
#sshfs hd_cd141@bwfor.cluster.uni-mannheim.de:/ bwFor  -o defer_permissions -o volname=Server
require(methods)
args = commandArgs(trailingOnly = TRUE)


version = "Transcript_info_v3"
StopFinder_version = "v12"
#on_desktop = "/Users/agkaessmann/bwFor"
on_desktop = ""
organism = args[1]
tissue = args[2]
core_num = as.numeric(args[3])
visualize = 0

#organism = "Human"
#tissue = "Brain"
#ORF_type = "pseudogene"

############################ THRESHOLDS #################################################

winsor_quantiles = c(0.1, 0.9)
ORF_length = 0
logit_threshold = 5
positive_codons_threshold = 4
covered_fraction_threshold = 0.1

#write(class(ORF_length),file = "/beegfs/work/ws/hd_cd141-RP-0/StopFinder/test")

########################### CLASSES ####################################################

setClass(
  "TranslationStop",
  slots = list(
    significant_stop = "numeric",
    significant_start = "numeric",
    LR = "numeric",
    frame =  "numeric",
    positive_codons =  "numeric",
    nonzero_codons =  "numeric",
    is_local_ORF = "numeric",
    covered_fraction =  "numeric",
    stop_entropy =  "numeric",
    with_ambigous = "numeric",
    coding_filtered = "numeric"
  )
)

setClass(
  "PredictionOutput",
  slots = list(
    significant_stop = "numeric",
    significant_start = "numeric",
    LR = "numeric",
    stop_entropy = "numeric",
    positive_codons = "numeric",
    frame = "numeric",
    covered_fraction = "numeric",
    with_ambigous = "numeric",
    is_local_ORF = "numeric",
    ORF_type = "character"
  )
)

############################ FUNCTIONS #################################################

winsorized_sum_uprolling = function(LR_codons_vector,
                                    ambigous_codons_vector,
                                    coding_codons_vector,
                                    take_ambigous) {
  if (length(LR_codons_vector)) {
    LR_codons_vector[is.na(LR_codons_vector)] = 0
    LR_uprolled = c()
    maxLR = 0
    maxLR_pos = 0
    full_maxLR = 0
    is_local_ORF = 0
    for (k in length(LR_codons_vector):1) {
      LR_codons_vector_cut = as.numeric(LR_codons_vector[k:length(LR_codons_vector)])
      ambigous_codons_vector_cut = ambigous_codons_vector[k:length(ambigous_codons_vector)]
      LR_codons_vector_cut_masked = c()
      if (length(coding_codons_vector) == 1) {
        LR_codons_vector_cut_masked = LR_codons_vector_cut[ambigous_codons_vector_cut ==
                                                             0]
      } else {
        coding_codons_vector_cut = coding_codons_vector[k:length(coding_codons_vector)]
        LR_codons_vector_cut_masked = LR_codons_vector_cut[ambigous_codons_vector_cut ==
                                                             0 & coding_codons_vector_cut == 0]
      }
      LR_full_sum = sum(LR_codons_vector_cut)
      if (LR_full_sum > full_maxLR) {
        full_maxLR = LR_full_sum
        maxLR_pos = k
      }
      lim = quantile(LR_codons_vector_cut_masked,
                     probs = winsor_quantiles,
                     type = 1)
      if (length(LR_codons_vector_cut_masked)) {
        replacement_low = 0
        replacement_high = 0
        if (length(which(LR_codons_vector_cut_masked > lim[1]))) {
          replacement_low = min(LR_codons_vector_cut_masked[which(LR_codons_vector_cut_masked >= lim[1])])
        }
        if (length(which(LR_codons_vector_cut_masked < lim[2]))) {
          replacement_high = max(LR_codons_vector_cut_masked[which(LR_codons_vector_cut_masked <= lim[2])])
        }
        LR_codons_vector_cut_masked[LR_codons_vector_cut_masked < lim[1]] = replacement_low
        LR_codons_vector_cut_masked[LR_codons_vector_cut_masked > lim[2]] = replacement_high
      }
      LR_sum = sum(LR_codons_vector_cut_masked)
      if (LR_sum > maxLR) {
        maxLR = LR_sum
        maxLR_pos = k
      }
      LR_uprolled = c(LR_uprolled, sum(LR_codons_vector_cut))
      if (ORF_length) {
        if (length(LR_codons_vector_cut_masked) >= ORF_length)
          break
      }
    }
    positive_codons_number = 0
    non_zero_codons_number = 0
    if (length(LR_uprolled)) {
      if (maxLR_pos == 0 & take_ambigous == "without_ambigous") {
        outer_LR = outer(rev(LR_uprolled), rev(LR_uprolled), FUN = "-")
        outer_LR[lower.tri(outer_LR)] = 0
        local_maxLR = max(outer_LR)
        start_end = which(outer_LR == local_maxLR, arr.ind = T)
        if (max(outer_LR) > 0 & start_end[2] - start_end[1] > 10) {
          maxLR = local_maxLR
          maxLR_pos = start_end[1]
          is_local_ORF = 1
        }
        positive_codons_number = sum(LR_codons_vector_cut_masked[maxLR_pos:start_end[2]] >
                                       0, na.rm = T)
        non_zero_codons_number = sum(LR_codons_vector_cut_masked[maxLR_pos:start_end[2]] !=
                                       0, na.rm = T)
      }
      else {
        positive_codons_number = sum(LR_codons_vector_cut_masked[maxLR_pos:length(LR_uprolled)] >
                                       0, na.rm = T)
        non_zero_codons_number = sum(LR_codons_vector_cut_masked[maxLR_pos:length(LR_uprolled)] !=
                                       0, na.rm = T)
      }
    }
    if (is.na(non_zero_codons_number)) {
      non_zero_codons_number = 0
    }
    if (is.na(positive_codons_number)) {
      positive_codons_number = 0
    }
    #try take maximal window
    return(
      c(
        maxLR,
        maxLR_pos,
        positive_codons_number,
        non_zero_codons_number,
        length(LR_codons_vector_cut_masked),
        is_local_ORF,
        rev(LR_uprolled)
      )
    )
  }
  else {
    return(c(0, 0, 0, 0, 0, 0))
  }
}


entropy = function(counts_codons_vector,
                   ambigous_codons_vector,
                   coding_codons_vector,
                   ORF_start) {
  counts_codons_vector_cut = counts_codons_vector[ORF_start:length(counts_codons_vector)]
  ambigous_codons_vector_cut = ambigous_codons_vector[ORF_start:length(ambigous_codons_vector)]
  counts_codons_vector_masked_cut = c()
  if (length(coding_codons_vector) == 1) {
    counts_codons_vector_masked_cut = counts_codons_vector_cut[ambigous_codons_vector_cut ==
                                                                 0]
  }
  else {
    coding_codons_vector_cut = coding_codons_vector[ORF_start:length(coding_codons_vector)]
    counts_codons_vector_masked_cut = counts_codons_vector_cut[ambigous_codons_vector_cut ==
                                                                 0 & coding_codons_vector_cut == 0]
  }
  if (length(counts_codons_vector_masked_cut) > 1) {
    if (sum(counts_codons_vector_masked_cut)) {
      codons_fractions = counts_codons_vector_masked_cut / sum(counts_codons_vector_masked_cut)
      return(sum(codons_fractions * log2(codons_fractions), na.rm = T) /
               log2(1 / length(counts_codons_vector_masked_cut)))
    }
    else {
      return(0)
    }
  }
  else {
    return("-1")
  }
}

split_tr_string = function(tr_string) {
  return(as.numeric(unlist(strsplit(
    as.character(tr_string), ","
  ))))
}

LR_by_codon = function(counts_matrix) {
  apply(counts_matrix, 2, function(x)
    (
      x[1] * log(3 * riboProb[1]) + x[2] * log(3 * riboProb[2]) + x[3] * log(3 * riboProb[3])
    ))
}

cut_mult3 = function(x) {
  return(x[1:(floor(length(x) / 3) * 3)])
}

sites2codons = function(sites_vector) {
  colSums(matrix(sites_vector, nrow = 3))
}

vector2matrix = function(counts_vector, frame) {
  matrix(cut_mult3(counts_vector[frame:length(counts_vector)]),
         nrow = 3,
         dimnames = list(c("f1", "f2", "f3"), 1:(length(
           cut_mult3(counts_vector[frame:length(counts_vector)])
         ) / 3)))
}

get_lwd_size = function(LR) {
  if (LR < 10) {
    return(2)
  }
  else if (LR < 30) {
    return(3)
  }
  else if (LR < 100) {
    return(4)
  }
  else if (LR < 1000) {
    return(5)
  }
  else {
    return(6)
  }
}

get_translated_ORFs = function(counts_vector,
                               coding_sites_vector,
                               take_ambigous,
                               take_coding) {
  new_stops = list()
  for (frame in c(1, 2, 3)) {
    new_stop = new("TranslationStop")
    stops_positions = c(frame,
                        which(stops_vector == frame),
                        length(stops_vector) - 3 + frame)
    if (length(stops_positions) > 1) {
      if (stops_positions[1] == stops_positions[2]) {
        stops_positions = stops_positions[-1]
      }
      if (length(stops_positions) > 1) {
        if (stops_positions[length(stops_positions)] == stops_positions[length(stops_positions) -
                                                                        1]) {
          stops_positions = stops_positions[1:(length(stops_positions) - 1)]
        }
      }
      if (length(stops_positions) > 1) {
        for (stop_num in 2:length(stops_positions)) {
          counts_vector_ORF = counts_vector[stops_positions[stop_num - 1]:(stops_positions[stop_num] -
                                                                             1)]
          counts_matrx_ORF = matrix(counts_vector_ORF, nrow = 3)
          LR = LR_by_codon(counts_matrx_ORF)
          ambigous_sites_vector_ORF = ambigous_sites_vector[stops_positions[stop_num -
                                                                              1]:(stops_positions[stop_num] - 1)]
          ambigous_sites_vector_ORF[is.na(ambigous_sites_vector_ORF)] = 0
          coding_sites_vector_ORF = coding_sites_vector[stops_positions[stop_num -
                                                                          1]:(stops_positions[stop_num] - 1)]
          ambigous_codons_vector_ORF = sites2codons(ambigous_sites_vector_ORF)
          coding_codons_vector_ORF = sites2codons(coding_sites_vector_ORF)
          if (take_ambigous == "with_ambigous") {
            ambigous_codons_vector_ORF = rep(0, length(ambigous_codons_vector_ORF))
            new_stop@with_ambigous = 1
          }
          else {
            new_stop@with_ambigous = 0
          }
          if (take_coding == "coding_not_filtered") {
            winsorized_sum_result = winsorized_sum_uprolling(LR,
                                                             ambigous_codons_vector_ORF,
                                                             "-1",
                                                             take_ambigous)
            new_stop@coding_filtered = 0
          }
          else{
            winsorized_sum_result = winsorized_sum_uprolling(
              LR,
              ambigous_codons_vector_ORF,
              coding_codons_vector_ORF,
              take_ambigous
            )
            new_stop@coding_filtered = 1
          }
          new_stop@LR = winsorized_sum_result[1]
          start_codon = winsorized_sum_result[2]
          new_stop@positive_codons = winsorized_sum_result[3]
          new_stop@nonzero_codons = winsorized_sum_result[4]
          all_codons = winsorized_sum_result[5]
          new_stop@is_local_ORF = winsorized_sum_result[6]
          new_stop@covered_fraction = 0
          if (all_codons) {
            new_stop@covered_fraction = new_stop@nonzero_codons / all_codons
          }
          if (new_stop@LR > 5 &
              new_stop@positive_codons > positive_codons_threshold &
              new_stop@covered_fraction > covered_fraction_threshold) {
            new_stop@significant_stop = stops_positions[stop_num] - 1
            new_stop@significant_start = stops_positions[stop_num - 1] + 3 +
              ((start_codon - 1) * 3)
            new_stop@frame = frame
            if (take_coding == "coding_not_filtered") {
              stop_entropy = entropy(counts_matrx_ORF,
                                     ambigous_codons_vector_ORF,
                                     "-1",
                                     start_codon)
            }
            else{
              stop_entropy = entropy(
                counts_matrx_ORF,
                ambigous_codons_vector_ORF,
                coding_codons_vector_ORF,
                start_codon
              )
            }
            new_stop@stop_entropy = stop_entropy
            new_stops = c(new_stops, new_stop)
          }
        }
      }
    }
  }
  return(new_stops)
}

predict_translation_four_times = function(counts_vector,coding_sites_vector) {
  #make first run without ambigous sites and coding-positions not filtered
  output_ORF_predictions = get_translated_ORFs(counts_vector,coding_sites_vector, "without_ambigous", "coding_not_filtered")

  #mask newly predicted ORFs
  if (length(output_ORF_predictions)) {
    for (ORF_num in (1:length(output_ORF_predictions))) {
      coding_sites_vector[output_ORF_predictions[[ORF_num]]@significant_start:output_ORF_predictions[[ORF_num]]@significant_stop] = 1
    }
  }

  #make 2nd run with coding positions here
  output_ORF_predictions = c(
    output_ORF_predictions,
    get_translated_ORFs(counts_vector,coding_sites_vector, "without_ambigous", "coding_filtered")
  )


  #handle ambigous positions
  #make 3rd run considering ambigous posistions
  amb_output_ORF_predictions = get_translated_ORFs(counts_vector, coding_sites_vector,"with_ambigous", "coding_not_filtered")
  output_ORF_predictions = c(output_ORF_predictions, amb_output_ORF_predictions)


  #filter newly predicted ambigous ORFs
  if (length(amb_output_ORF_predictions)) {
    for (ORF_num in (1:length(amb_output_ORF_predictions))) {
      coding_sites_vector[amb_output_ORF_predictions[[ORF_num]]@significant_start:amb_output_ORF_predictions[[ORF_num]]@significant_stop] = 1
    }
  }

  #make 4th run considering ambigous posistions and newly predicted coding posistions
  output_ORF_predictions = c(
    output_ORF_predictions,
    get_translated_ORFs(counts_vector,coding_sites_vector, "with_ambigous", "coding_filtered")
  )
  return(output_ORF_predictions)
}

prepare_prediction_output = function(ORF_predictions,transcript_type) {
  new_output = new("PredictionOutput")
  if (length(ORF_predictions)) {
    for (significant_stop_num in (1:length(ORF_predictions))) {
      if (sum(ORF_predictions[[significant_stop_num]]@significant_stop == new_output@significant_stop) == 0) {
        for (output_slot in c(
          "significant_stop",
          "significant_start",
          "LR",
          "stop_entropy",
          "positive_codons",
          "frame",
          "covered_fraction",
          "with_ambigous",
          "is_local_ORF"
        )) {
          slot(new_output, output_slot) = c(slot(new_output, output_slot),
                                            slot(ORF_predictions[[significant_stop_num]], output_slot))
        }
      }
      else {
        if ((new_output@significant_start)[which(new_output@significant_stop == ORF_predictions[[significant_stop_num]]@significant_stop)] >
            ORF_predictions[[significant_stop_num]]@significant_start) {
          slot(new_output, "significant_start")[which(new_output@significant_stop == ORF_predictions[[significant_stop_num]]@significant_stop)] =
            ORF_predictions[[significant_stop_num]]@significant_start
        }
      }
    }

    #make ORF classification!!!!!!!
    new_output@ORF_type = rep("x", length(new_output@significant_stop))
    if (transcript_type == "protein_coding" &
        length(new_output@significant_stop)) {
      canonical_stop = transcripts[line_num, "canonical_stop"] - 1
      canonical_start = transcripts[line_num, "canonical_start"]
      if (is.na(canonical_stop)) {
        canonical_stop = (new_output@significant_stop)[which.max(new_output@significant_stop - new_output@significant_start)]
      }
      stop_num = which(new_output@significant_stop == canonical_stop)
      slot(new_output, "ORF_type")[stop_num] = "c"
      #take position of the previous stop
      stops_positions = c((new_output@frame)[stop_num],
                          which(stops_vector == (new_output@frame)[stop_num]),
                          length(stops_vector) - 3 + (new_output@frame)[stop_num]
      )
      previous_stop = stops_positions[which(stops_positions == canonical_stop) - 1]
      slot(new_output, "ORF_type")[which(
        new_output@significant_stop > canonical_stop &
          new_output@significant_start > canonical_stop
      )] = "dn"
      slot(new_output, "ORF_type")[which(
        new_output@significant_stop > canonical_stop &
          new_output@significant_start < canonical_stop
      )] = "do"
      slot(new_output, "ORF_type")[which(
        new_output@significant_stop < canonical_start &
          new_output@significant_start < canonical_start
      )] = "un"
      slot(new_output, "ORF_type")[which(
        new_output@significant_stop < canonical_stop &
          new_output@significant_stop > canonical_start &
          new_output@significant_start < canonical_start
      )] = "uo"
    }
  }

  #check coding overlap
  check_ORF_vector = c()
  if (length(new_output@ORF_type)) {
    new_coding_sites_vector = coding_sites_vector
    for (ORF_num in 1:length(new_output@ORF_type)) {
      if ((new_output@ORF_type)[ORF_num] == "c") {
        ORF_range = (new_output@significant_start)[ORF_num]:(new_output@significant_stop)[ORF_num]
        new_coding_sites_vector[matrix(cut_mult3(ORF_range),nrow=3)[1,]] = 1
      }
    }
    for (ORF_num in 1:length(new_output@ORF_type)) {
      if ((new_output@ORF_type)[ORF_num] == "x") {
        ORF_range = (new_output@significant_start)[ORF_num]:(new_output@significant_stop)[ORF_num]
        if (sum(new_coding_sites_vector[matrix(cut_mult3(ORF_range),nrow=3)[1,]])) {
          check_ORF_vector = c(check_ORF_vector, F)
        }
        else{
          check_ORF_vector = c(check_ORF_vector, T)
        }
      }
      else if ((new_output@ORF_type)[ORF_num] == "c") {
        check_ORF_vector = c(check_ORF_vector, T)
      }
      else {
        ORF_range = (new_output@significant_start)[ORF_num]:(new_output@significant_stop)[ORF_num]###!!!! check here
        if((new_output@ORF_type)[ORF_num] == "un" || (new_output@ORF_type)[ORF_num] == "uo") {
          if((new_output@significant_start)[ORF_num] > min(which(new_coding_sites_vector>0))) {
            check_ORF_vector = c(check_ORF_vector, F)
          }
          else {
            check_ORF_vector = c(check_ORF_vector, T)
          }
        }
        if((new_output@ORF_type)[ORF_num] == "dn" || (new_output@ORF_type)[ORF_num] == "do") {
          if((new_output@significant_stop)[ORF_num] < max(which(new_coding_sites_vector>0))) {
            check_ORF_vector = c(check_ORF_vector, F)
          }
          else {
            check_ORF_vector = c(check_ORF_vector, T)
          }
        }
      }
    }
  }


  for (output_slot in c(
    "significant_stop",
    "significant_start",
    "LR",
    "stop_entropy",
    "positive_codons",
    "frame",
    "covered_fraction",
    "with_ambigous",
    "is_local_ORF",
    "ORF_type"
  )) {
    slot(new_output, output_slot) = slot(new_output, output_slot)[check_ORF_vector] # filter non-canonical ORFs overlap
  }
  return(new_output)
}

#ADD SIMON's POISSON!!!

############################ MAIN BODY #################################################
#get file with fraction estimates
allRiboProbs = read.table(
  paste(
    on_desktop,
    "/beegfs/work/ws/hd_cd141-RP-0/StopFinder/allRiboProbs",
    sep = ""
  ),
  sep = "\t"
)

#uniform probs and ribo probs
unifProb = c(1 / 3, 1 / 3, 1 / 3)
riboProb = as.numeric(allRiboProbs[allRiboProbs$organism == organism &
                                     allRiboProbs$tissue == tissue, c("f1", "f2", "f3")])


transcripts = read.table(
  paste(
    on_desktop,
    paste(
      "/beegfs/work/ws/hd_cd141-RP-0/StopFinder",
      version,
      organism,
      tissue,
      sep = "/"
    ),
    "/",
    "all_transcripts_add_",
    core_num,
    ".txt",
    sep = ""
  ),
  header = TRUE,
  sep = "\t",
  stringsAsFactors = FALSE
)



for (line_num in 1:nrow(transcripts)) {
  #for (line_num in 1:100) {
  #ambigous sites: 1 - ambigous, 0 - unique
  ambigous_sites_vector = split_tr_string(transcripts[line_num, "ambigous_sites"])
  #overlap with coding sites
  coding_sites_vector = split_tr_string(transcripts[line_num, "coding_sites"])
  tr_biotype = transcripts[line_num,"tr_biotype"]
  if(is.na(tr_biotype)) {
    tr_biotype = "uknown"
  }
  else if(tr_biotype == "") {
    tr_biotype = "uknown"
  }
  ###take read counts: 1 - overlap with coding, 0 - no overlap with coding
  #RP
  RP_counts_vector = split_tr_string(transcripts[line_num, "RP_counts"])
  #turn to matrix
  RP_counts_matrix_f1 = vector2matrix(RP_counts_vector, 1)
  RP_counts_matrix_f2 = vector2matrix(RP_counts_vector, 2)
  RP_counts_matrix_f3 = vector2matrix(RP_counts_vector, 3)

  #TR
  TR_counts_vector = split_tr_string(transcripts[line_num, "TR_counts"])
  #turn to matrix
  TR_counts_matrix_f1 = vector2matrix(TR_counts_vector, 1)
  TR_counts_matrix_f2 = vector2matrix(TR_counts_vector, 2)
  TR_counts_matrix_f3 = vector2matrix(TR_counts_vector, 3)

  ###locate stop-codons
  RNA_seq = as.character(transcripts[line_num, "RNA_sequence"])
  if (length(RNA_seq) == 0) next
  stops_vector_linear = cut_mult3(c(rep(0, nchar(RNA_seq))))
  stops_vector_linear[unlist(gregexpr("TGA|TAG|TAA", RNA_seq))] = 1
  stops_vector = stops_vector_linear * c(1, 2, 3)

  #calculate maximal winsorized sum of likelihood ratios by rolling upwards
  RP_ORF_predictions = predict_translation_four_times(RP_counts_vector,coding_sites_vector)
  TR_ORF_predictions = predict_translation_four_times(TR_counts_vector,coding_sites_vector)

  ###calculate whole-transcript LR
  #RP
  RP_LR_f1 = LR_by_codon(RP_counts_matrix_f1)
  RP_LR_f2 = LR_by_codon(RP_counts_matrix_f2)
  RP_LR_f3 = LR_by_codon(RP_counts_matrix_f3)
  RP_LR_matrix = matrix(c(RP_LR_f1[1:length(RP_LR_f3)], RP_LR_f2[1:length(RP_LR_f3)], RP_LR_f3),
                        nrow = 3,
                        byrow = T)

  #TR
  TR_LR_f1 = LR_by_codon(TR_counts_matrix_f1)
  TR_LR_f2 = LR_by_codon(TR_counts_matrix_f2)
  TR_LR_f3 = LR_by_codon(TR_counts_matrix_f3)
  TR_LR_matrix = matrix(c(TR_LR_f1[1:length(TR_LR_f3)], TR_LR_f2[1:length(TR_LR_f3)], TR_LR_f3),
                        nrow = 3,
                        byrow = T)


  #newly predicted coding
  predicted_coding = matrix(NA, nrow = 3, ncol = length(coding_sites_vector))
  if (length(RP_ORF_predictions)) {
    for (significant_stop_num in (1:length(RP_ORF_predictions))) {
      predicted_coding[RP_ORF_predictions[[significant_stop_num]]@frame,
                       RP_ORF_predictions[[significant_stop_num]]@significant_start:(RP_ORF_predictions[[significant_stop_num]]@significant_stop -
                                                                                       3)] =  1
    }
  }
  predicted_coding = predicted_coding * c(3, 2, 1)


  #predicted as coding for TR (false positives)
  TR_predicted_coding = matrix(NA, nrow = 3, ncol = length(coding_sites_vector))
  if (length(TR_ORF_predictions)) {
    for (significant_stop_num in (1:length(TR_ORF_predictions))) {
      TR_predicted_coding[TR_ORF_predictions[[significant_stop_num]]@frame,
                          TR_ORF_predictions[[significant_stop_num]]@significant_start:(TR_ORF_predictions[[significant_stop_num]]@significant_stop -
                                                                                          3)] =  1
    }
  }
  TR_predicted_coding = TR_predicted_coding * c(3, 2, 1)

  RP_prediction_output = prepare_prediction_output(RP_ORF_predictions,tr_biotype)
  TR_prediction_output = prepare_prediction_output(TR_ORF_predictions,tr_biotype)

  #prepare the output
  for (output_category in c(
    "ORF_type",
    "significant_stop",
    "significant_start",
    "frame",
    "LR",
    "with_ambigous",
    "covered_fraction",
    "stop_entropy",
    "positive_codons",
    "is_local_ORF"
  )) {
    transcripts[line_num, paste("RP_prediction", output_category, sep = "_")] = paste(slot(RP_prediction_output, output_category), collapse = ",")
    transcripts[line_num, paste("TR_prediction", output_category, sep = "_")] = paste(slot(TR_prediction_output, output_category), collapse = ",")
  }

  #visualize
  if (visualize) {
    par(mfrow = c(6, 1), mar = c(3.1, 3.1, 3.1, 3.1))
    #visualize RP
    barplot(
      RP_counts_matrix_f1,
      col = c(rep(
        c("red", "green", "blue"), ncol(RP_counts_matrix_f1)
      )),
      border = NA,
      beside = T,
      ylim = c(0, quantile(RP_counts_matrix_f1, probs = 0.95))
    )
    barplot(
      RP_LR_matrix,
      col = c(rep(
        c("red", "green", "blue"), ncol(RP_counts_matrix_f1)
      )),
      border = NA,
      beside = T,
      ylim = c(-20, 20)
    )
    stops_vector_linear[stops_vector_linear == 0] = NA
    plot(
      stops_vector_linear * c(3, 2, 1),
      pch = "|",
      col = c("red", "green", "blue"),
      cex = 3,
      ylim = c(0, 4),
      bty = "n",
      yaxt = "n"
    )
    colors_vector = c("red", "green", "blue")
    if (length(RP_ORF_predictions)) {
      for (significant_stop_num in (1:length(RP_ORF_predictions))) {
        to_plot = rep(NA, length(predicted_coding[RP_ORF_predictions[[significant_stop_num]]@frame, ]))
        range_to_plot = RP_ORF_predictions[[significant_stop_num]]@significant_start:(RP_ORF_predictions[[significant_stop_num]]@significant_stop -
                                                                                        3)
        to_plot[range_to_plot] = predicted_coding[RP_ORF_predictions[[significant_stop_num]]@frame, range_to_plot]
        lwd_size = get_lwd_size(RP_ORF_predictions[[significant_stop_num]]@LR)
        lty_value = 1
        if (RP_ORF_predictions[[significant_stop_num]]@with_ambigous) {
          lty_value = 3
        }
        lines(to_plot,
              col = colors_vector[RP_ORF_predictions[[significant_stop_num]]@frame],
              lwd = lwd_size,
              lty = lty_value)
      }
    }

    #visualize TR
    #par(mfrow=c(3,1),mar=c(3.1,3.1,3.1,3.1))

    barplot(
      TR_counts_matrix_f1,
      col = c(rep(
        c("red", "green", "blue"), ncol(TR_counts_matrix_f1)
      )),
      border = NA,
      beside = T,
      ylim = c(0, quantile(TR_counts_matrix_f1, probs = 0.95))
    )
    barplot(
      TR_LR_matrix,
      col = c(rep(
        c("red", "green", "blue"), ncol(TR_counts_matrix_f1)
      )),
      border = NA,
      beside = T,
      ylim = c(-20, 20)
    )
    stops_vector_linear[stops_vector_linear == 0] = NA
    plot(
      stops_vector_linear * c(3, 2, 1),
      pch = "|",
      col = c("red", "green", "blue"),
      cex = 3,
      ylim = c(0, 4),
      bty = "n",
      yaxt = "n"
    )
    colors_vector = c("red", "green", "blue")
    if (length(TR_ORF_predictions)) {
      for (significant_stop_num in (1:length(TR_ORF_predictions))) {
        to_plot = rep(NA, length(predicted_coding[TR_ORF_predictions[[significant_stop_num]]@frame, ]))
        range_to_plot = TR_ORF_predictions[[significant_stop_num]]@significant_start:(TR_ORF_predictions[[significant_stop_num]]@significant_stop -
                                                                                        3)
        to_plot[range_to_plot] = TR_predicted_coding[TR_ORF_predictions[[significant_stop_num]]@frame, range_to_plot]
        lwd_size = get_lwd_size(TR_ORF_predictions[[significant_stop_num]]@LR)
        lty_value = 1
        if (TR_ORF_predictions[[significant_stop_num]]@with_ambigous) {
          lty_value = 3
        }
        lines(to_plot,
              col = colors_vector[TR_ORF_predictions[[significant_stop_num]]@frame],
              lwd = lwd_size,
              lty = lty_value)
      }
    }
  }
}

print_col_names = T

write.table(
  transcripts,
  file = paste(
    on_desktop,
    paste(
      "/beegfs/work/ws/hd_cd141-RP-0/StopFinder",
      version,
      organism,
      tissue,
      sep = "/"
    ),
    "/",
    "stop_predictions_",
    StopFinder_version,
    "_RP_TR_all_",
    core_num,
    ".txt",
#    print_suffix,
    sep = ""
  ),
  quote = F,
  sep = "\t",
  row.names = F,
  col.names = print_col_names
)
