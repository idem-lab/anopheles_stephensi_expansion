# code from David Duncan
library(rgbif)

gbif_data  <- occ_download(
  pred('taxonKey', 1),
  pred_in('basisOfRecord', 
          c("MACHINE_OBSERVATION", "HUMAN_OBSERVATION")),
  pred('country', "AU"),
  #pred_within('POLYGON((112.5 -44.5, 154.5 -44.5, 154.5 -9.5, 112.5 -9.5, 112.5 -44.5))'),
  pred('hasGeospatialIssue', "FALSE"),
  pred('occurrenceStatus', "PRESENT"),
  pred_lt("coordinateUncertaintyInMeters",1000),
  pred_gte('year', 2000),
  pred_or(pred_not(
    pred_in("establishmentMeans",
            c("MANAGED","INTRODUCED"))),
    pred_isnull("establishmentMeans")),
  format = "SIMPLE_CSV")

gbif_animalia <- occ_download_get('0430467-210914110416597') %>%
  occ_download_import() %>% 
  write_csv(
    glue('{data_path}gbif_ala_occurrence/gbif_animalia_{today()}.csv')
  )