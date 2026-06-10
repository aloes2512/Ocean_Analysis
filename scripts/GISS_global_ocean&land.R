global(ocean_mask, fun = "sum", na.rm = TRUE) / ncell(ocean_mask)  # Should be ~60-70%
plot(ocean_mask, main = "Ocean Mask Coverage")
#====
compareGeom(gistemp, ocean_mask, stopOnError = FALSE)  # Should show identical=TRUE
res(gistemp); res(ocean_mask)  # Resolutions must match
#======
# First time slice only
test_slice <- gistemp[[1]]
masked_slice <- mask(test_slice, ocean_mask)
#global(masked_slice, "mean", na.rm = TRUE)  # Should differ from unmasked
#global(test_slice,"mean",na.rm=T)
# more debug steps
# Check actual values in your mask
values(ocean_mask)[1:20]  # Should see many TRUE/FALSE
table(values(ocean_mask))  # TRUE proportion ~60-70% for oceans
#global(ocean_mask, "sum") / ncell(ocean_mask)  # Coverage fraction
# fixing 1
#gistemp_ocean <- mask(gistemp, ocean_mask, inverse = TRUE)  # Keeps where mask=FALSE (land)
# OR for ocean-only:
#ocean_slice <- mask(test_slice, ocean_mask)  # Should now work ✓
#global(ocean_slice,"mean",na.rm=T)
# fix 2
land_mask <- ocean_mask == FALSE  # TRUE = land cells to mask OUT
ocean_slice <- mask(test_slice, land_mask)
global(land_mask ,"mean",na.rm=T    )
global(!land_mask,"mean",na.rm=T)
# apply
# Define from land_mask first (global(land_mask, "mean") = 0.38 = 38% land coverage ✓)
ocean_mask <- !land_mask                    # 62% ocean coverage ✓
all_equal(ocean_mask,!land_mask)
# Apply masks correctly:
gistemp_land  <- mask(gistemp, land_mask)   # Keep where land_mask=TRUE (land cells)
gistemp_ocean <- mask(gistemp, ocean_mask)  # Keep where ocean_mask=TRUE (ocean cells)

# Spatial mean for each time layer (returns data.frame)
ts_ocean <- global(gistemp_ocean, fun = "mean", na.rm = TRUE)
ts_land<-global(gistemp_land, fun = "mean", na.rm = TRUE)
range(ts_ocean)
# Extract temperature values and dates
tibble(dates = time(gistemp_ocean), # POSIX dates from NetCDF,
ocean_anom = ts_ocean$mean ) %>% # Vector of global ocean means
ggplot(aes(x=dates,y=ocean_anom))+geom_line()
range(ts_ocean) #-1.147346  1.998164
range(ts_land)
# avoid logical completely
# Convert to numeric 0/1 (GUARANTEED to work)
land_mask_num  <- as.numeric(land_mask)   # 1 = land, 0 = ocean
ocean_mask_num <- 1 - land_mask_num       # 0 = land, 1 = ocean ✓
summary(land_mask_num)
summary(ocean_mask_num)
table(values(land_mask_num) > 0.9)
table(values(ocean_mask_num)>0.9)
land_mask_final <- land_mask_num >= 0.9   # Catches 1.0 exactly
ocean_mask_final <- ocean_mask_num >= 0.9
# Verify arithmetic works
global(land_mask_final, "sum") / ncell(land_mask_num)   # ~0.38
global(ocean_mask_final, "sum") / ncell(land_mask_num) # ~0.62

# Mask with numeric masks
gistemp_land  <- mask(gistemp, land_mask_final)   # Keep land
gistemp_ocean <- mask(gistemp, ocean_mask_final)  # Keep ocean
global(gistemp_land[[100]], "mean", na.rm=TRUE)    # Land warmer
global(gistemp_ocean[[100]], "mean", na.rm=TRUE)   # Ocean cooler

ts_land  <- global(gistemp_land, "mean", na.rm=TRUE)
ts_ocean <- global(gistemp_ocean, "mean", na.rm=TRUE)

range(ts_ocean$mean)   # ±0.15°C ✓
range(ts_land$mean)    # ±0.30°C ✓

crs(gistemp, describe=TRUE)
crs(land_mask, describe=TRUE)

compareGeom(gistemp, land_mask, stopOnError=FALSE, crs=TRUE)
# Must return TRUE for crs, extent, AND resolution

gistemp_land1  <- mask(gistemp[[1]], land_mask_final)   # Keep land
gistemp_ocean1 <- mask(gistemp[[1]], ocean_mask_final)  # Keep ocean

global(gistemp_land1,"mean",na.rm=T)
global(gistemp_ocean1,"mean",na.rm=T)
