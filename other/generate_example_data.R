# Need to ensure your file and object names match
load('data/amph_afr_df.rda')
amph_afr_df <- input
save(amph_afr_df, file = 'data/amph_afr_df.rda', compress = 'xz')

load('data/kruger_park.rda')
kruger_park <- target_area
save(kruger_park, file = 'data/kruger_park.rda')  # no compression, destroys obj
