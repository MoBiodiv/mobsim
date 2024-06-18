# snapshot test with Thoms process

    Code
      jitter_species(simdat, seed = 42L)
    Output
      $census
                  x         y    species
      1  0.70562431 0.7492636 species_01
      2  0.66944759 0.4828926 species_01
      3  0.93533619 0.4812669 species_02
      4  0.96695333 0.1083010 species_02
      5  0.91843733 0.6411184 species_03
      6  0.94340404 0.5321446 species_07
      7  0.06576721 0.4925216 species_08
      8  0.43193198 1.0121692 species_08
      9  0.28507829 0.7226997 species_09
      10 0.84556285 0.1318787 species_10
      
      $x_min_max
      [1] 0 1
      
      $y_min_max
      [1] 0 1
      
      attr(,"class")
      [1] "community"

# snapshot test with Poisson process

    Code
      jitter_species(simdat, seed = 42L)
    Output
      $census
                 x         y    species
      1  0.9285156 0.4567952 species_01
      2  0.9507850 0.7181657 species_01
      3  0.2804926 0.9548565 species_02
      4  0.8248006 0.2756131 species_02
      5  0.6453768 0.4616657 species_03
      6  0.5254246 0.9530632 species_07
      7  0.7406310 1.0010929 species_08
      8  0.1387093 0.1403538 species_08
      9  0.6559310 0.4611085 species_09
      10 0.7201800 0.5575449 species_10
      
      $x_min_max
      [1] 0 1
      
      $y_min_max
      [1] 0 1
      
      attr(,"class")
      [1] "community"

# snapshot test wide jitter

    Code
      jitter_species(simdat, sd = 3, seed = 42L)
    Output
      $census
                   x          y    species
      1   4.80479007  0.4662330 species_01
      2   4.76861335  0.1998621 species_01
      3  -0.75311134  6.5163538 species_02
      4  -0.72149420  6.1433879 species_02
      5   2.00419128  0.4536032 species_03
      6   2.83566323  4.4337049 species_07
      7   1.27452950  7.3295913 species_08
      8   1.64069426  7.8492389 species_08
      9  -0.03223401 -3.4299938 species_09
      10  5.36501362 -0.7016997 species_10
      
      $x_min_max
      [1] 0 1
      
      $y_min_max
      [1] 0 1
      
      attr(,"class")
      [1] "community"

# snapshot test with singletons only

    Code
      jitter_species(simdat, seed = 42L)
    Output
      $census
                 x          y    species
      1  0.9285156 0.47079047 species_01
      2  0.9314284 0.74197871 species_02
      3  0.2897708 0.92078364 species_03
      4  0.8367763 0.25264094 species_04
      5  0.6457882 0.46095961 species_05
      6  0.5180347 0.94637403 species_06
      7  0.7517035 0.97538390 species_07
      8  0.1337200 0.09092281 species_08
      9  0.6771765 0.45059241 species_09
      10 0.7044376 0.57353388 species_10
      
      $x_min_max
      [1] 0 1
      
      $y_min_max
      [1] 0 1
      
      attr(,"class")
      [1] "community"

