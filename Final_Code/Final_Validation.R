### Load in the clean reference table and the clean BEAST/BFAST table
### Reference size threshold (static parameter) 
### Create an accuracy table
  ### For-Loop, for each sample
      ### For-Loop, for each year:
        ### True Positive: if (reference size > reference threshold and reference size and beast size both positive) OR (reference size > reference thd for following year and beast size positive)
        ### True Negative : if reference size =< reference threshold and beast size = 0
        ### False Positive: beast size positive, reference size =< reference threshold
        ### False Negative: if beast size = 0, reference size > reference threshold
        ### NA (No Data)

