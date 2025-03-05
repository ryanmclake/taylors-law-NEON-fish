# Gererate a pass taxa Count table to estimate abundance

# clear your memory
gc()

# set an exclusion function
'%!in%' <- function(x,y)!('%in%'(x,y))

# identify lake sites to remove
lakesites <- c("PRPO", "PRLA", "CRAM", "LIRO", "TOOK")

# create a passTaxaCount table to quantify 3-pass depletion abundance
bulk_for_passTaxaCount <- fsh_bulkCount %>% select(eventID, boutEndDate, passNumber, taxonID, namedLocation, bulkFishCount) %>%
  group_by(eventID, boutEndDate, passNumber, namedLocation, taxonID) %>%
  summarize(bulk_sum = sum(bulkFishCount))

fixed_reaches <- fsh_fieldData %>%
  select(namedLocation, siteID, fixedRandomReach, domainID) %>% 
  group_by(namedLocation) %>%
  slice(1)

pass_taxa_count_table <- fsh_perFish %>% select(eventID, boutEndDate, passNumber, namedLocation, taxonID) %>%
  group_by(eventID, boutEndDate, passNumber, namedLocation, taxonID) %>%
  tally() %>%
  rename(fish_sum = n) %>%
  left_join(., bulk_for_passTaxaCount, by = c("eventID", "boutEndDate", "passNumber", "taxonID", "namedLocation")) |>
  mutate(bulk_sum = ifelse(is.na(bulk_sum), 0, bulk_sum),
         taxa_total = bulk_sum + fish_sum) |>
  ungroup()|>
  select(-bulk_sum, -fish_sum) %>%
  left_join(., fixed_reaches, by = "namedLocation") %>%
  filter(fixedRandomReach == "fixed") %>%
  filter(siteID %!in% lakesites)

joined_fish <- pass_taxa_count_table %>%
  dplyr::group_by(eventID, boutEndDate, domainID, siteID, namedLocation, taxonID) %>%
  pivot_wider(names_from = passNumber, values_from = taxa_total, values_fill = 0) %>%
  ungroup() %>%
  mutate(flag = ifelse(`2`>`1` & `3`>`1`, "1", "0")) %>%
  mutate(flag2 = ifelse(`2`== 0 & `3`>0, "1", "0")) %>%
  mutate(flag3 = ifelse(`2`== 0 & `3`== 0, "1", "0")) %>%
  filter(flag == "0") %>%
  filter(flag2 == "0") %>%
  filter(flag3 == "0") %>%
  select(-flag, -flag2, -flag3) %>%
  mutate(row = row_number(.)) %>%
  write_csv(., "output/fsh_passTaxaCount.csv")
