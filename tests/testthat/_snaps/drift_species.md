# snapshot test with drift

    Code
      drift_x_species(simdat, drift = 0.1)
    Output
      $census
                 x         y    species
      1  1.0148060 0.4577418 species_01
      2  1.0370754 0.7191123 species_02
      3  0.3861395 0.9346722 species_03
      4  0.9304476 0.2554288 species_04
      5  0.7417455 0.4622928 species_05
      6  0.6190959 0.9400145 species_06
      7  0.8365883 0.9782264 species_07
      8  0.2346666 0.1174874 species_08
      9  0.7569923 0.4749971 species_09
      10 0.8050648 0.5603327 species_10
      
      $x_min_max
      [1] 0 1
      
      $y_min_max
      [1] 0 1
      
      attr(,"class")
      [1] "community"

---

    Code
      drift_y_species(simdat, drift = 0.1)
    Output
      $census
                 x         y    species
      1  0.9148060 0.5577418 species_01
      2  0.9370754 0.8191123 species_02
      3  0.2861395 1.0346722 species_03
      4  0.8304476 0.3554288 species_04
      5  0.6417455 0.5622928 species_05
      6  0.5190959 1.0400145 species_06
      7  0.7365883 1.0782264 species_07
      8  0.1346666 0.2174874 species_08
      9  0.6569923 0.5749971 species_09
      10 0.7050648 0.6603327 species_10
      
      $x_min_max
      [1] 0 1
      
      $y_min_max
      [1] 0 1
      
      attr(,"class")
      [1] "community"

