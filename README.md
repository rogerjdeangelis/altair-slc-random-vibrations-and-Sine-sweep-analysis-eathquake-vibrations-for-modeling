# altair-slc-random-vibrations-and-Sine-sweep-analysis-eathquake-vibrations-for-modeling
Altair slc random vibrations and Sine sweep analysis eathquake vibrations for modeling
    %let pgm=altair-slc-random-vibrations-and-Sine-sweep-analysis-eathquake-vibrations-for-modeling;

    %stop_submission;

    Altair slc random vibrations and Sine sweep analysis eathquake vibrations for modeling

    Too long to post here,see
    https://github.com/rogerjdeangelis/altair-slc-random-vibrations-and-Sine-sweep-analysis-eathquake-vibrations-for-modeling

    Plots
    https://github.com/rogerjdeangelis/altair-slc-random-vibrations-and-Sine-sweep-analysis-eathquake-vibrations-for-modeling/blob/main/random_vibration_psd.png
    https://github.com/rogerjdeangelis/altair-slc-random-vibrations-and-Sine-sweep-analysis-eathquake-vibrations-for-modeling/blob/main/random_vibration_time_history.png
    https://github.com/rogerjdeangelis/altair-slc-random-vibrations-and-Sine-sweep-analysis-eathquake-vibrations-for-modeling/blob/main/sine_sweep_frequency.png


    This analysis is essential for ensuring products survive real-world vibration environments
    from
     your smartphone surviving a drop,
     to a car handling a bumpy road,
     to a rocket staying intact during launch,
     to a building surviving a eathquake

    The code generates the vibration profiles that engineers use to test
    and validate products before they reach consumers.

    This is the profile of the random vibrations

    Random Vibrations

    Frequency    PSD       Interpretation
    (Hz)        (g²/Hz)

    20          0.01        Low energy at very low frequencies
    80          0.03        Energy increases 3× (more vibration!)
    200         0.03        Energy stays high (sustained vibration)
    500         0.005       Energy drops 6× (high frequencies damped)

    /*___
    |  _ \   _ __  _ __ ___   __ _ _ __ __ _ _ __ ___
    | |_) | | `_ \| `__/ _ \ / _` | `__/ _` | `_ ` _ \
    |  _ <  | |_) | | | (_) | (_| | | | (_| | | | | | |
    |_| \_\ | .__/|_|  \___/ \__, |_|  \__,_|_| |_| |_|
            |_|              |___/
    */

    options validvarname=v7;
    options set=RHOME "C:\Progra~1\R\R-4.5.2\bin\r";
    proc r;
    submit;
    # ============================================
    # RANDOM VIBRATION AND SINE SWEEP ANALYSIS
    # Converted from Python to R - CORRECTED VERSION
    # ============================================

    # Load required libraries
    library(tidyverse)
    library(signal)
    library(pracma)
    library(ggplot2)
    library(gridExtra)

    # -----------------------------
    # Set working directory
    # -----------------------------
    setwd("d:/wpswrkx")

    # Create directories for output if they don't exist
    if (!dir.exists("d:/csv")) dir.create("d:/csv", recursive = TRUE)
    if (!dir.exists("d:/png")) dir.create("d:/png", recursive = TRUE)

    # -----------------------------
    # Random vibration example
    # -----------------------------
    fs <- 2000
    T <- 20.0
    dt <- 1 / fs
    t <- seq(0, T - dt, by = dt)
    n <- length(t)

    # PSD breakpoint profile
    f_bp <- c(20, 80, 200, 500)
    psd_bp <- c(0.01, 0.03, 0.03, 0.005)

    # Print parameters
    cat("========================================\n")
    cat("RANDOM VIBRATION ANALYSIS\n")
    cat("========================================\n")
    cat("Frequency breakpoints (Hz):", f_bp, "\n")
    cat("PSD values (g²/Hz):", psd_bp, "\n")
    cat("Sampling frequency:", fs, "Hz\n")
    cat("Duration:", T, "seconds\n")
    cat("Number of samples:", n, "\n")

    # -----------------------------
    # Overall GRMS from area under PSD curve
    # -----------------------------
    psd_area <- trapz(f_bp, psd_bp)
    target_grms <- sqrt(psd_area)

    cat("PSD area:", psd_area, "g²\n")
    cat("Target GRMS:", target_grms, "g\n")

    # -----------------------------
    # Create a synthetic random signal
    # -----------------------------
    set.seed(7)
    random_vibration <- rnorm(n, mean = 0, sd = 1)

    cat("\nFirst 10 random values:\n")
    cat(head(random_vibration, 10), "\n")

    # -----------------------------
    # Apply low-pass filter
    # -----------------------------
    b <- c(1.0)
    a <- c(1.0, -0.98)
    filtered <- filter(b, a, random_vibration)

    # Scale to match target GRMS
    current_std <- sd(filtered)
    random_vibration_scaled <- filtered / current_std * target_grms

    cat("\nRandom vibration statistics after filtering and scaling:\n")
    cat("Mean:", mean(random_vibration_scaled), "\n")
    cat("Std Dev:", sd(random_vibration_scaled), "\n")
    cat("Min:", min(random_vibration_scaled), "\n")
    cat("Max:", max(random_vibration_scaled), "\n")

    # -----------------------------
    # Estimate PSD using Welch's method
    # -----------------------------
    # Use the specgram function from signal package
    nperseg <- 2048

    # Method 1: Using specgram (works in most signal package versions)
    tryCatch({
        # Calculate spectrogram with Welch averaging
        spec_result <- specgram(random_vibration_scaled,
                                n = nperseg,
                                Fs = fs,
                                overlap = nperseg / 2)

        # Average the spectrogram to get PSD (Welch-like)
        # The spectrogram gives power in each time-frequency bin
        # Average across time to get the PSD estimate
        pxx <- rowMeans(Mod(spec_result$S)^2)
        f_est <- spec_result$f

        # Trim to positive frequencies (usually specgram gives up to fs/2)
        # The PSD from specgram is already scaled appropriately
        cat("\nPSD estimation using specgram completed.\n")
        cat("Frequency range:", min(f_est), "to", max(f_est), "Hz\n")
        cat("Number of frequency bins:", length(f_est), "\n")

    }, error = function(e) {
        # Method 2: Fallback - simple periodogram if specgram fails
        cat("\nWarning: specgram failed, using simple periodogram.\n")
        cat("Error message:", e$message, "\n")

        # Use FFT directly for a simple periodogram
        nfft <- 2048
        fft_result <- fft(random_vibration_scaled[1:nfft])
        pxx <- (Mod(fft_result)[1:(nfft/2 + 1)])^2 / (fs * nfft)
        f_est <- (0:(nfft/2)) * fs / nfft
    })

    # If specgram failed, pxx and f_est might not exist
    if (!exists("pxx") || !exists("f_est")) {
        cat("\nUsing alternative PSD estimation method...\n")

        # Use pspectrum if available (newer signal package)
        if (exists("pspectrum")) {
            psd_result <- pspectrum(random_vibration_scaled, fs = fs,
                                    nfft = nperseg, window = hanning(nperseg))
            f_est <- psd_result$freq
            pxx <- psd_result$spec
        } else {
            # Manual Welch method
            nseg <- floor(n / nperseg)
            pxx_accum <- rep(0, nperseg/2 + 1)

            for (seg in 1:nseg) {
                start_idx <- (seg - 1) * nperseg + 1
                end_idx <- start_idx + nperseg - 1
                segment <- random_vibration_scaled[start_idx:end_idx]

                # Apply Hanning window
                window <- 0.5 * (1 - cos(2 * pi * (0:(nperseg-1)) / (nperseg - 1)))
                segment_windowed <- segment * window

                # FFT
                fft_result <- fft(segment_windowed)
                pxx_seg <- (Mod(fft_result)[1:(nperseg/2 + 1)])^2 / (fs * sum(window^2))
                pxx_accum <- pxx_accum + pxx_seg
            }

            pxx <- pxx_accum / nseg
            f_est <- (0:(nperseg/2)) * fs / nperseg
            cat("Manual Welch PSD estimation completed.\n")
        }
    }

    # Calculate estimated GRMS
    estimated_grms <- sqrt(trapz(f_est, pxx))

    cat("\nEstimated GRMS from generated signal:", estimated_grms, "g\n")
    cat("Difference from target:", abs(estimated_grms - target_grms), "g\n")

    # -----------------------------
    # Sine sweep example
    # -----------------------------
    cat("\n========================================\n")
    cat("SINE SWEEP ANALYSIS\n")
    cat("========================================\n")

    f0 <- 10
    f1 <- 300
    k <- log(f1 / f0) / T

    phase <- 2 * pi * f0 * (exp(k * t) - 1) / k
    sine_sweep <- sin(phase)
    inst_freq <- f0 * exp(k * t)

    cat("Start frequency:", f0, "Hz\n")
    cat("End frequency:", f1, "Hz\n")
    cat("Sweep rate:", k, "\n")
    cat("First 10 time points:\n")
    cat(head(t, 10), "\n")

    # -----------------------------
    # Create data frames for export
    # -----------------------------
    example_timeseries <- data.frame(
        time_s = t,
        random_vibration_g = random_vibration_scaled,
        sine_sweep = sine_sweep,
        inst_freq_hz = inst_freq
    )

    random_psd_profile <- data.frame(
        freq_hz = f_bp,
        psd_g2_per_hz = psd_bp
    )

    random_psd_estimate <- data.frame(
        freq_hz = f_est,
        psd_est_g2_per_hz = pxx
    )

    # -----------------------------
    # Save data to CSV
    # -----------------------------
    write.csv(example_timeseries,
              file = "d:/csv/example_timeseries.csv",
              row.names = FALSE)

    write.csv(random_psd_profile,
              file = "d:/csv/random_psd_profile.csv",
              row.names = FALSE)

    write.csv(random_psd_estimate,
              file = "d:/csv/random_psd_estimate.csv",
              row.names = FALSE)

    cat("\nData saved to CSV files in d:/csv/\n")

    # -----------------------------
    # Create plots using ggplot2
    # -----------------------------

    # Plot 1: Random Vibration PSD
    plot1 <- ggplot() +
        geom_point(data = random_psd_profile,
                   aes(x = freq_hz, y = psd_g2_per_hz, color = "Target PSD"),
                   size = 3) +
        geom_line(data = random_psd_profile,
                  aes(x = freq_hz, y = psd_g2_per_hz, color = "Target PSD"),
                  linewidth = 1) +
        geom_line(data = random_psd_estimate,
                  aes(x = freq_hz, y = psd_est_g2_per_hz, color = "Welch estimate"),
                  alpha = 0.8, linewidth = 1) +
        scale_x_log10(breaks = c(1, 10, 100, 1000),
                      labels = c("1", "10", "100", "1000")) +
        scale_y_log10() +
        labs(x = "Frequency (Hz)",
             y = "PSD (g²/Hz)",
             title = paste("Random Vibration PSD",
                           "target GRMS =", round(target_grms, 3),
                           "estimated GRMS =", round(estimated_grms, 3)),
             color = "") +
        theme_minimal() +
        theme(legend.position = c(0.85, 0.85),
              panel.grid.minor = element_line(linetype = "dotted")) +
        scale_color_manual(values = c("Target PSD" = "blue",
                                      "Welch estimate" = "red"))

    # Save plot
    ggsave("d:/png/random_vibration_psd.png",
           plot = plot1, width = 8, height = 4.5, dpi = 160)

    # Plot 2: Sine Sweep Signal
    plot2 <- ggplot(example_timeseries, aes(x = time_s, y = sine_sweep)) +
        geom_line(color = "red", linewidth = 0.5) +
        labs(x = "Time (s)",
             y = "Amplitude",
             title = "Exponential Sine Sweep") +
        theme_minimal() +
        theme(panel.grid.minor = element_line(linetype = "dotted"))

    ggsave("d:/png/sine_sweep_signal.png",
           plot = plot2, width = 8, height = 4.5, dpi = 160)

    # Plot 3: Instantaneous Frequency
    plot3 <- ggplot(example_timeseries, aes(x = time_s, y = inst_freq_hz)) +
        geom_line(color = "darkred", linewidth = 1) +
        labs(x = "Time (s)",
             y = "Frequency (Hz)",
             title = "Sine Sweep Instantaneous Frequency") +
        theme_minimal() +
        theme(panel.grid.minor = element_line(linetype = "dotted"))

    ggsave("d:/png/sine_sweep_frequency.png",
           plot = plot3, width = 8, height = 4.5, dpi = 160)

    # Plot 4: Random Vibration Time History (subset)
    plot4 <- ggplot(example_timeseries[1:1000, ],
                    aes(x = time_s, y = random_vibration_g)) +
        geom_line(color = "blue", linewidth = 0.5) +
        labs(x = "Time (s)",
             y = "Acceleration (g)",
             title = "Random Vibration Time History (First 1000 samples)") +
        theme_minimal() +
        theme(panel.grid.minor = element_line(linetype = "dotted"))

    ggsave("d:/png/random_vibration_time_history.png",
           plot = plot4, width = 8, height = 4.5, dpi = 160)

    # -----------------------------
    # Print summary statistics
    # -----------------------------
    cat("\n========================================\n")
    cat("SUMMARY STATISTICS\n")
    cat("========================================\n")

    cat("\nRandom Vibration Statistics:\n")
    cat("Mean:", mean(random_vibration_scaled), "\n")
    cat("Std Dev:", sd(random_vibration_scaled), "\n")
    cat("Min:", min(random_vibration_scaled), "\n")
    cat("Max:", max(random_vibration_scaled), "\n")
    cat("GRMS (target):", target_grms, "\n")
    cat("GRMS (estimated):", estimated_grms, "\n")

    cat("\nSine Sweep Statistics:\n")
    cat("Mean:", mean(sine_sweep), "\n")
    cat("Std Dev:", sd(sine_sweep), "\n")
    cat("Min:", min(sine_sweep), "\n")
    cat("Max:", max(sine_sweep), "\n")
    cat("Frequency range:", min(inst_freq), "to", max(inst_freq), "Hz\n")

    cat("\n========================================\n")
    cat("Analysis Complete!\n")
    cat("========================================\n")
    endsubmit;
    run;

    /*           _               _
      ___  _   _| |_ _ __  _   _| |_
     / _ \| | | | __| `_ \| | | | __|
    | (_) | |_| | |_| |_) | |_| | |_
     \___/ \__,_|\__| .__/ \__,_|\__|
                    |_|
    */

    /**************************************************************************************************************************/
    /*  PLOTS                                                                                                                 */
    /*                                                                                                                        */
    /* random_vibration_time_history.png                                                                                      */
    /* random_vibration_psd.png                                                                                               */
    /* sine_sweep_frequency.png                                                                                               */
    /* random_vibration_psd.png                                                                                               */
    /*                                                                                                                        */
    /* Altair SLC                                                                                                             */
    /*                                                                                                                        */
    /* ========================================                                                                               */
    /* RANDOM VIBRATION ANALYSIS                                                                                              */
    /* ========================================                                                                               */
    /* Frequency breakpoints (Hz): 20 80 200 500                                                                              */
    /* PSD values (gÂ²/Hz): 0.01 0.03 0.03 0.005                                                                              */
    /* Sampling frequency: 2000 Hz                                                                                            */
    /* Duration: 20 seconds                                                                                                   */
    /* Number of samples: 40000                                                                                               */
    /* PSD area: 10.05 gÂ²                                                                                                    */
    /* Target GRMS: 3.170173 g                                                                                                */
    /* First 10 random values:                                                                                                */
    /* 2.287247 -1.196772 -0.6942925 -0.412293 -0.9706733 -0.9472799 0.7481393 -0.1169552 0.1526576 2.189978                  */
    /* Random vibration statistics after filtering and scaling:                                                               */
    /* Mean: 0.09765129                                                                                                       */
    /* Std Dev: 3.170173                                                                                                      */
    /* Min: -11.07822                                                                                                         */
    /* Max: 10.42558                                                                                                          */
    /* PSD estimation using specgram completed.                                                                               */
    /* Frequency range: 0 to 999.0234 Hz                                                                                      */
    /* Number of frequency bins: 1024                                                                                         */
    /* Estimated GRMS from generated signal: 2780.572 g                                                                       */
    /* Difference from target: 2777.402 g                                                                                     */
    /* ========================================                                                                               */
    /* SINE SWEEP ANALYSIS                                                                                                    */
    /* ========================================                                                                               */
    /* Start frequency: 10 Hz                                                                                                 */
    /* End frequency: 300 Hz                                                                                                  */
    /* Sweep rate: 0.1700599                                                                                                  */
    /* First 10 time points:                                                                                                  */
    /* 0 5e-04 0.001 0.0015 0.002 0.0025 0.003 0.0035 0.004 0.0045                                                            */
    /* Data saved to CSV files in d:/csv/                                                                                     */
    /* ========================================                                                                               */
    /* SUMMARY STATISTICS                                                                                                     */
    /* ========================================                                                                               */
    /* Random Vibration Statistics:                                                                                           */
    /* Mean: 0.09765129                                                                                                       */
    /* Std Dev: 3.170173                                                                                                      */
    /* Min: -11.07822                                                                                                         */
    /* Max: 10.42558                                                                                                          */
    /* GRMS (target): 3.170173                                                                                                */
    /* GRMS (estimated): 2780.572                                                                                             */
    /* Sine Sweep Statistics:                                                                                                 */
    /* Mean: 0.0007883135                                                                                                     */
    /* Std Dev: 0.7071077                                                                                                     */
    /* Min: -1                                                                                                                */
    /* Max: 1                                                                                                                 */
    /* Frequency range: 10 to 299.9745 Hz                                                                                     */
    /* ========================================                                                                               */
    /* Analysis Complete!                                                                                                     */
    /* ========================================                                                                               */
    /**************************************************************************************************************************/

    PLOTS

    random_vibration_time_history.png
    random_vibration_psd.png
    sine_sweep_frequency.png
    random_vibration_psd.png

    Altair SLC

    ========================================
    RANDOM VIBRATION ANALYSIS
    ========================================
    Frequency breakpoints (Hz): 20 80 200 500
    PSD values (gÂ²/Hz): 0.01 0.03 0.03 0.005
    Sampling frequency: 2000 Hz
    Duration: 20 seconds
    Number of samples: 40000
    PSD area: 10.05 gÂ²
    Target GRMS: 3.170173 g
    First 10 random values:
    2.287247 -1.196772 -0.6942925 -0.412293 -0.9706733 -0.9472799 0.7481393 -0.1169552 0.1526576 2.189978
    Random vibration statistics after filtering and scaling:
    Mean: 0.09765129
    Std Dev: 3.170173
    Min: -11.07822
    Max: 10.42558
    PSD estimation using specgram completed.
    Frequency range: 0 to 999.0234 Hz
    Number of frequency bins: 1024
    Estimated GRMS from generated signal: 2780.572 g
    Difference from target: 2777.402 g
    ========================================
    SINE SWEEP ANALYSIS
    ========================================
    Start frequency: 10 Hz
    End frequency: 300 Hz
    Sweep rate: 0.1700599
    First 10 time points:
    0 5e-04 0.001 0.0015 0.002 0.0025 0.003 0.0035 0.004 0.0045
    Data saved to CSV files in d:/csv/
    ========================================
    SUMMARY STATISTICS
    ========================================
    Random Vibration Statistics:
    Mean: 0.09765129
    Std Dev: 3.170173
    Min: -11.07822
    Max: 10.42558
    GRMS (target): 3.170173
    GRMS (estimated): 2780.572
    Sine Sweep Statistics:
    Mean: 0.0007883135
    Std Dev: 0.7071077
    Min: -1
    Max: 1
    Frequency range: 10 to 299.9745 Hz
    ========================================
    Analysis Complete!
    ========================================

    /*
    | | ___   __ _
    | |/ _ \ / _` |
    | | (_) | (_| |
    |_|\___/ \__, |
             |___/
    */

    1                                          Altair SLC         16:03 Thursday, June 25, 2026

    NOTE: Copyright 2002-2025 World Programming, an Altair Company
    NOTE: Altair SLC 2026 (05.26.01.00.000758)
          Licensed to Roger DeAngelis
    NOTE: This session is executing on the X64_WIN11PRO platform and is running in 64 bit mode

    NOTE: AUTOEXEC processing beginning; file is C:\wpsoto\autoexec.sas
    NOTE: Library workx assigned as follows:
          Engine:        SAS7BDAT
          Physical Name: d:\wpswrkx

    NOTE: Library wpdx assigned as follows:
          Engine:        WPD
          Physical Name: d:\wpswrkx

    NOTE: Library slchelp assigned as follows:
          Engine:        WPD
          Physical Name: C:\Progra~1\Altair\SLC\2026\sashelp


    LOG:  16:03:57
    NOTE: 1 record was written to file PRINT

    NOTE: The data step took :
          real time : 0.047
          cpu time  : 0.031


    NOTE: Format num2mis output
    NOTE: Format $chr2mis output
    NOTE: Procedure format step took :
          real time : 0.015
          cpu time  : 0.000


    NOTE: AUTOEXEC processing completed

    1          options validvarname=v7;
    2         options set=RHOME "C:\Progra~1\R\R-4.5.2\bin\r";
    3         proc r;
    4         submit;
    5         # ============================================
    6         # RANDOM VIBRATION AND SINE SWEEP ANALYSIS
    7         # Converted from Python to R - CORRECTED VERSION
    8         # ============================================
    9
    10        # Load required libraries
    11        library(tidyverse)
    12        library(signal)
    13        library(pracma)
    14        library(ggplot2)
    15        library(gridExtra)
    16
    17        # -----------------------------
    18        # Set working directory
    19        # -----------------------------
    20        setwd("d:/wpswrkx")
    21
    22        # Create directories for output if they don't exist
    23        if (!dir.exists("d:/csv")) dir.create("d:/csv", recursive = TRUE)
    24        if (!dir.exists("d:/png")) dir.create("d:/png", recursive = TRUE)
    25
    26        # -----------------------------
    27        # Random vibration example
    28        # -----------------------------
    29        fs <- 2000
    30        T <- 20.0
    31        dt <- 1 / fs
    32        t <- seq(0, T - dt, by = dt)
    33        n <- length(t)
    34
    35        # PSD breakpoint profile
    36        f_bp <- c(20, 80, 200, 500)
    37        psd_bp <- c(0.01, 0.03, 0.03, 0.005)
    38
    39        # Print parameters
    40        cat("========================================\n")
    41        cat("RANDOM VIBRATION ANALYSIS\n")
    42        cat("========================================\n")
    43        cat("Frequency breakpoints (Hz):", f_bp, "\n")
    44        cat("PSD values (gÂ²/Hz):", psd_bp, "\n")
    45        cat("Sampling frequency:", fs, "Hz\n")
    46        cat("Duration:", T, "seconds\n")
    47        cat("Number of samples:", n, "\n")
    48
    49        # -----------------------------
    50        # Overall GRMS from area under PSD curve
    51        # -----------------------------
    52        psd_area <- trapz(f_bp, psd_bp)
    53        target_grms <- sqrt(psd_area)
    54
    55        cat("PSD area:", psd_area, "gÂ²\n")
    56        cat("Target GRMS:", target_grms, "g\n")
    57
    58        # -----------------------------
    59        # Create a synthetic random signal
    60        # -----------------------------
    61        set.seed(7)
    62        random_vibration <- rnorm(n, mean = 0, sd = 1)
    63
    64        cat("\nFirst 10 random values:\n")
    65        cat(head(random_vibration, 10), "\n")
    66
    67        # -----------------------------
    68        # Apply low-pass filter
    69        # -----------------------------
    70        b <- c(1.0)
    71        a <- c(1.0, -0.98)
    72        filtered <- filter(b, a, random_vibration)
    73
    74        # Scale to match target GRMS
    75        current_std <- sd(filtered)
    76        random_vibration_scaled <- filtered / current_std * target_grms
    77
    78        cat("\nRandom vibration statistics after filtering and scaling:\n")
    79        cat("Mean:", mean(random_vibration_scaled), "\n")
    80        cat("Std Dev:", sd(random_vibration_scaled), "\n")
    81        cat("Min:", min(random_vibration_scaled), "\n")
    82        cat("Max:", max(random_vibration_scaled), "\n")
    83
    84        # -----------------------------
    85        # Estimate PSD using Welch's method
    86        # -----------------------------
    87        # Use the specgram function from signal package
    88        nperseg <- 2048
    89
    90        # Method 1: Using specgram (works in most signal package versions)
    91        tryCatch({
    92            # Calculate spectrogram with Welch averaging
    93            spec_result <- specgram(random_vibration_scaled,
    94                                    n = nperseg,
    95                                    Fs = fs,
    96                                    overlap = nperseg / 2)
    97
    98            # Average the spectrogram to get PSD (Welch-like)
    99            # The spectrogram gives power in each time-frequency bin
    100           # Average across time to get the PSD estimate
    101           pxx <- rowMeans(Mod(spec_result$S)^2)
    102           f_est <- spec_result$f
    103
    104           # Trim to positive frequencies (usually specgram gives up to fs/2)
    105           # The PSD from specgram is already scaled appropriately
    106           cat("\nPSD estimation using specgram completed.\n")
    107           cat("Frequency range:", min(f_est), "to", max(f_est), "Hz\n")
    108           cat("Number of frequency bins:", length(f_est), "\n")
    109
    110       }, error = function(e) {
    111           # Method 2: Fallback - simple periodogram if specgram fails
    112           cat("\nWarning: specgram failed, using simple periodogram.\n")
    113           cat("Error message:", e$message, "\n")
    114
    115           # Use FFT directly for a simple periodogram
    116           nfft <- 2048
    117           fft_result <- fft(random_vibration_scaled[1:nfft])
    118           pxx <- (Mod(fft_result)[1:(nfft/2 + 1)])^2 / (fs * nfft)
    119           f_est <- (0:(nfft/2)) * fs / nfft
    120       })
    121
    122       # If specgram failed, pxx and f_est might not exist
    123       if (!exists("pxx") || !exists("f_est")) {
    124           cat("\nUsing alternative PSD estimation method...\n")
    125
    126           # Use pspectrum if available (newer signal package)
    127           if (exists("pspectrum")) {
    128               psd_result <- pspectrum(random_vibration_scaled, fs = fs,
    129                                       nfft = nperseg, window = hanning(nperseg))
    130               f_est <- psd_result$freq
    131               pxx <- psd_result$spec
    132           } else {
    133               # Manual Welch method
    134               nseg <- floor(n / nperseg)
    135               pxx_accum <- rep(0, nperseg/2 + 1)
    136
    137               for (seg in 1:nseg) {
    138                   start_idx <- (seg - 1) * nperseg + 1
    139                   end_idx <- start_idx + nperseg - 1
    140                   segment <- random_vibration_scaled[start_idx:end_idx]
    141
    142                   # Apply Hanning window
    143                   window <- 0.5 * (1 - cos(2 * pi * (0:(nperseg-1)) / (nperseg - 1)))
    144                   segment_windowed <- segment * window
    145
    146                   # FFT
    147                   fft_result <- fft(segment_windowed)
    148                   pxx_seg <- (Mod(fft_result)[1:(nperseg/2 + 1)])^2 / (fs * sum(window^2))
    149                   pxx_accum <- pxx_accum + pxx_seg
    150               }
    151
    152               pxx <- pxx_accum / nseg
    153               f_est <- (0:(nperseg/2)) * fs / nperseg
    154               cat("Manual Welch PSD estimation completed.\n")
    155           }
    156       }
    157
    158       # Calculate estimated GRMS
    159       estimated_grms <- sqrt(trapz(f_est, pxx))
    160
    161       cat("\nEstimated GRMS from generated signal:", estimated_grms, "g\n")
    162       cat("Difference from target:", abs(estimated_grms - target_grms), "g\n")
    163
    164       # -----------------------------
    165       # Sine sweep example
    166       # -----------------------------
    167       cat("\n========================================\n")
    168       cat("SINE SWEEP ANALYSIS\n")
    169       cat("========================================\n")
    170
    171       f0 <- 10
    172       f1 <- 300
    173       k <- log(f1 / f0) / T
    174
    175       phase <- 2 * pi * f0 * (exp(k * t) - 1) / k
    176       sine_sweep <- sin(phase)
    177       inst_freq <- f0 * exp(k * t)
    178
    179       cat("Start frequency:", f0, "Hz\n")
    180       cat("End frequency:", f1, "Hz\n")
    181       cat("Sweep rate:", k, "\n")
    182       cat("First 10 time points:\n")
    183       cat(head(t, 10), "\n")
    184
    185       # -----------------------------
    186       # Create data frames for export
    187       # -----------------------------
    188       example_timeseries <- data.frame(
    189           time_s = t,
    190           random_vibration_g = random_vibration_scaled,
    191           sine_sweep = sine_sweep,
    192           inst_freq_hz = inst_freq
    193       )
    194
    195       random_psd_profile <- data.frame(
    196           freq_hz = f_bp,
    197           psd_g2_per_hz = psd_bp
    198       )
    199
    200       random_psd_estimate <- data.frame(
    201           freq_hz = f_est,
    202           psd_est_g2_per_hz = pxx
    203       )
    204
    205       # -----------------------------
    206       # Save data to CSV
    207       # -----------------------------
    208       write.csv(example_timeseries,
    209                 file = "d:/csv/example_timeseries.csv",
    210                 row.names = FALSE)
    211
    212       write.csv(random_psd_profile,
    213                 file = "d:/csv/random_psd_profile.csv",
    214                 row.names = FALSE)
    215
    216       write.csv(random_psd_estimate,
    217                 file = "d:/csv/random_psd_estimate.csv",
    218                 row.names = FALSE)
    219
    220       cat("\nData saved to CSV files in d:/csv/\n")
    221
    222       # -----------------------------
    223       # Create plots using ggplot2
    224       # -----------------------------
    225
    226       # Plot 1: Random Vibration PSD
    227       plot1 <- ggplot() +
    228           geom_point(data = random_psd_profile,
    229                      aes(x = freq_hz, y = psd_g2_per_hz, color = "Target PSD"),
    230                      size = 3) +
    231           geom_line(data = random_psd_profile,
    232                     aes(x = freq_hz, y = psd_g2_per_hz, color = "Target PSD"),
    233                     linewidth = 1) +
    234           geom_line(data = random_psd_estimate,
    235                     aes(x = freq_hz, y = psd_est_g2_per_hz, color = "Welch estimate"),
    236                     alpha = 0.8, linewidth = 1) +
    237           scale_x_log10(breaks = c(1, 10, 100, 1000),
    238                         labels = c("1", "10", "100", "1000")) +
    239           scale_y_log10() +
    240           labs(x = "Frequency (Hz)",
    241                y = "PSD (gÂ²/Hz)",
    242                title = paste("Random Vibration PSD",
    243                              "target GRMS =", round(target_grms, 3),
    244                              "estimated GRMS =", round(estimated_grms, 3)),
    245                color = "") +
    246           theme_minimal() +
    247           theme(legend.position = c(0.85, 0.85),
    248                 panel.grid.minor = element_line(linetype = "dotted")) +
    249           scale_color_manual(values = c("Target PSD" = "blue",
    250                                         "Welch estimate" = "red"))
    251
    252       # Save plot
    253       ggsave("d:/png/random_vibration_psd.png",
    254              plot = plot1, width = 8, height = 4.5, dpi = 160)
    255
    256       # Plot 2: Sine Sweep Signal
    257       plot2 <- ggplot(example_timeseries, aes(x = time_s, y = sine_sweep)) +
    258           geom_line(color = "red", linewidth = 0.5) +
    259           labs(x = "Time (s)",
    260                y = "Amplitude",
    261                title = "Exponential Sine Sweep") +
    262           theme_minimal() +
    263           theme(panel.grid.minor = element_line(linetype = "dotted"))
    264
    265       ggsave("d:/png/sine_sweep_signal.png",
    266              plot = plot2, width = 8, height = 4.5, dpi = 160)
    267
    268       # Plot 3: Instantaneous Frequency
    269       plot3 <- ggplot(example_timeseries, aes(x = time_s, y = inst_freq_hz)) +
    270           geom_line(color = "darkred", linewidth = 1) +
    271           labs(x = "Time (s)",
    272                y = "Frequency (Hz)",
    273                title = "Sine Sweep Instantaneous Frequency") +
    274           theme_minimal() +
    275           theme(panel.grid.minor = element_line(linetype = "dotted"))
    276
    277       ggsave("d:/png/sine_sweep_frequency.png",
    278              plot = plot3, width = 8, height = 4.5, dpi = 160)
    279
    280       # Plot 4: Random Vibration Time History (subset)
    281       plot4 <- ggplot(example_timeseries[1:1000, ],
    282                       aes(x = time_s, y = random_vibration_g)) +
    283           geom_line(color = "blue", linewidth = 0.5) +
    284           labs(x = "Time (s)",
    285                y = "Acceleration (g)",
    286                title = "Random Vibration Time History (First 1000 samples)") +
    287           theme_minimal() +
    288           theme(panel.grid.minor = element_line(linetype = "dotted"))
    289
    290       ggsave("d:/png/random_vibration_time_history.png",
    291              plot = plot4, width = 8, height = 4.5, dpi = 160)
    292
    293       # -----------------------------
    294       # Print summary statistics
    295       # -----------------------------
    296       cat("\n========================================\n")
    297       cat("SUMMARY STATISTICS\n")
    298       cat("========================================\n")
    299
    300       cat("\nRandom Vibration Statistics:\n")
    301       cat("Mean:", mean(random_vibration_scaled), "\n")
    302       cat("Std Dev:", sd(random_vibration_scaled), "\n")
    303       cat("Min:", min(random_vibration_scaled), "\n")
    304       cat("Max:", max(random_vibration_scaled), "\n")
    305       cat("GRMS (target):", target_grms, "\n")
    306       cat("GRMS (estimated):", estimated_grms, "\n")
    307
    308       cat("\nSine Sweep Statistics:\n")
    309       cat("Mean:", mean(sine_sweep), "\n")
    310       cat("Std Dev:", sd(sine_sweep), "\n")
    311       cat("Min:", min(sine_sweep), "\n")
    312       cat("Max:", max(sine_sweep), "\n")
    313       cat("Frequency range:", min(inst_freq), "to", max(inst_freq), "Hz\n")
    314
    315       cat("\n========================================\n")
    316       cat("Analysis Complete!\n")
    317       cat("========================================\n")
    318       endsubmit;
    NOTE: Using R version 4.5.2 (2025-10-31 ucrt) from C:\Program Files\R\R-4.5.2

    NOTE: Submitting statements to R:

    > # ============================================
    > # RANDOM VIBRATION AND SINE SWEEP ANALYSIS
    > # Converted from Python to R - CORRECTED VERSION
    > # ============================================
    >
    > # Load required libraries
    > library(tidyverse)
    -- Attaching core tidyverse packages ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- tidyverse 2.0.0 --
    v dplyr     1.2.1     v readr     2.2.0
    v forcats   1.0.1     v stringr   1.6.0
    v ggplot2   4.0.2     v tibble    3.3.1
    v lubridate 1.9.5     v tidyr     1.3.2
    v purrr     1.2.1
    -- Conflicts ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- tidyverse_conflicts() --
    x dplyr::filter() masks stats::filter()
    x dplyr::lag()    masks stats::lag()
    i Use the conflicted package (<http://conflicted.r-lib.org/>) to force all conflicts to become errors
    Warning messages:
    1: package 'ggplot2' was built under R version 4.5.3
    2: package 'tidyr' was built under R version 4.5.3
    3: package 'dplyr' was built under R version 4.5.3
    4: package 'stringr' was built under R version 4.5.3
    > library(signal)
    Attaching package: 'signal'
    The following object is masked from 'package:dplyr':
        filter
    The following objects are masked from 'package:stats':
        filter, poly
    Warning message:
    package 'signal' was built under R version 4.5.3
    > library(pracma)
    Attaching package: 'pracma'
    The following objects are masked from 'package:signal':
        conv, ifft, interp1, pchip, polyval, roots
    The following object is masked from 'package:purrr':
        cross
    Warning message:
    package 'pracma' was built under R version 4.5.3
    > library(ggplot2)
    > library(gridExtra)
    Attaching package: 'gridExtra'
    The following object is masked from 'package:dplyr':
        combine
    >
    > # -----------------------------
    > # Set working directory
    > # -----------------------------
    > setwd("d:/wpswrkx")
    >
    > # Create directories for output if they don't exist
    > if (!dir.exists("d:/csv")) dir.create("d:/csv", recursive = TRUE)
    > if (!dir.exists("d:/png")) dir.create("d:/png", recursive = TRUE)
    >
    > # -----------------------------
    > # Random vibration example
    > # -----------------------------
    > fs <- 2000
    > T <- 20.0
    > dt <- 1 / fs
    > t <- seq(0, T - dt, by = dt)
    > n <- length(t)
    >
    > # PSD breakpoint profile
    > f_bp <- c(20, 80, 200, 500)
    > psd_bp <- c(0.01, 0.03, 0.03, 0.005)
    >
    > # Print parameters
    > cat("========================================\n")
    > cat("RANDOM VIBRATION ANALYSIS\n")
    > cat("========================================\n")
    > cat("Frequency breakpoints (Hz):", f_bp, "\n")
    > cat("PSD values (gÂ²/Hz):", psd_bp, "\n")
    > cat("Sampling frequency:", fs, "Hz\n")
    > cat("Duration:", T, "seconds\n")
    > cat("Number of samples:", n, "\n")
    >
    > # -----------------------------
    > # Overall GRMS from area under PSD curve
    > # -----------------------------
    > psd_area <- trapz(f_bp, psd_bp)
    > target_grms <- sqrt(psd_area)
    >
    > cat("PSD area:", psd_area, "gÂ²\n")
    > cat("Target GRMS:", target_grms, "g\n")
    >
    > # -----------------------------
    > # Create a synthetic random signal
    > # -----------------------------
    > set.seed(7)
    > random_vibration <- rnorm(n, mean = 0, sd = 1)
    >
    > cat("\nFirst 10 random values:\n")
    > cat(head(random_vibration, 10), "\n")
    >
    > # -----------------------------
    > # Apply low-pass filter
    > # -----------------------------
    > b <- c(1.0)
    > a <- c(1.0, -0.98)
    > filtered <- filter(b, a, random_vibration)
    >
    > # Scale to match target GRMS
    > current_std <- sd(filtered)
    > random_vibration_scaled <- filtered / current_std * target_grms
    >
    > cat("\nRandom vibration statistics after filtering and scaling:\n")
    > cat("Mean:", mean(random_vibration_scaled), "\n")
    > cat("Std Dev:", sd(random_vibration_scaled), "\n")
    > cat("Min:", min(random_vibration_scaled), "\n")
    > cat("Max:", max(random_vibration_scaled), "\n")
    >
    > # -----------------------------
    > # Estimate PSD using Welch's method
    > # -----------------------------
    > # Use the specgram function from signal package
    > nperseg <- 2048
    >
    > # Method 1: Using specgram (works in most signal package versions)
    > tryCatch({
    +     # Calculate spectrogram with Welch averaging
    +     spec_result <- specgram(random_vibration_scaled,
    +                             n = nperseg,
    +                             Fs = fs,
    +                             overlap = nperseg / 2)
    +
    +     # Average the spectrogram to get PSD (Welch-like)
    +     # The spectrogram gives power in each time-frequency bin
    +     # Average across time to get the PSD estimate
    +     pxx <- rowMeans(Mod(spec_result$S)^2)
    +     f_est <- spec_result$f
    +
    +     # Trim to positive frequencies (usually specgram gives up to fs/2)
    +     # The PSD from specgram is already scaled appropriately
    +     cat("\nPSD estimation using specgram completed.\n")
    +     cat("Frequency range:", min(f_est), "to", max(f_est), "Hz\n")
    +     cat("Number of frequency bins:", length(f_est), "\n")
    +
    + }, error = function(e) {
    +     # Method 2: Fallback - simple periodogram if specgram fails
    +     cat("\nWarning: specgram failed, using simple periodogram.\n")

    2                                                                                                                         Altair SLC

    +     cat("Error message:", e$message, "\n")
    +
    +     # Use FFT directly for a simple periodogram
    +     nfft <- 2048
    +     fft_result <- fft(random_vibration_scaled[1:nfft])
    +     pxx <- (Mod(fft_result)[1:(nfft/2 + 1)])^2 / (fs * nfft)
    +     f_est <- (0:(nfft/2)) * fs / nfft
    + })
    >
    > # If specgram failed, pxx and f_est might not exist
    > if (!exists("pxx") || !exists("f_est")) {
    +     cat("\nUsing alternative PSD estimation method...\n")
    +
    +     # Use pspectrum if available (newer signal package)
    +     if (exists("pspectrum")) {
    +         psd_result <- pspectrum(random_vibration_scaled, fs = fs,
    +                                 nfft = nperseg, window = hanning(nperseg))
    +         f_est <- psd_result$freq
    +         pxx <- psd_result$spec
    +     } else {
    +         # Manual Welch method
    +         nseg <- floor(n / nperseg)
    +         pxx_accum <- rep(0, nperseg/2 + 1)
    +
    +         for (seg in 1:nseg) {
    +             start_idx <- (seg - 1) * nperseg + 1
    +             end_idx <- start_idx + nperseg - 1
    +             segment <- random_vibration_scaled[start_idx:end_idx]
    +
    +             # Apply Hanning window
    +             window <- 0.5 * (1 - cos(2 * pi * (0:(nperseg-1)) / (nperseg - 1)))
    +             segment_windowed <- segment * window
    +
    +             # FFT
    +             fft_result <- fft(segment_windowed)
    +             pxx_seg <- (Mod(fft_result)[1:(nperseg/2 + 1)])^2 / (fs * sum(window^2))
    +             pxx_accum <- pxx_accum + pxx_seg
    +         }
    +
    +         pxx <- pxx_accum / nseg
    +         f_est <- (0:(nperseg/2)) * fs / nperseg
    +         cat("Manual Welch PSD estimation completed.\n")
    +     }
    + }
    >
    > # Calculate estimated GRMS
    > estimated_grms <- sqrt(trapz(f_est, pxx))
    >
    > cat("\nEstimated GRMS from generated signal:", estimated_grms, "g\n")
    > cat("Difference from target:", abs(estimated_grms - target_grms), "g\n")
    >
    > # -----------------------------
    > # Sine sweep example
    > # -----------------------------
    > cat("\n========================================\n")
    > cat("SINE SWEEP ANALYSIS\n")
    > cat("========================================\n")
    >
    > f0 <- 10
    > f1 <- 300
    > k <- log(f1 / f0) / T
    >
    > phase <- 2 * pi * f0 * (exp(k * t) - 1) / k
    > sine_sweep <- sin(phase)
    > inst_freq <- f0 * exp(k * t)
    >
    > cat("Start frequency:", f0, "Hz\n")
    > cat("End frequency:", f1, "Hz\n")
    > cat("Sweep rate:", k, "\n")
    > cat("First 10 time points:\n")
    > cat(head(t, 10), "\n")
    >
    > # -----------------------------
    > # Create data frames for export
    > # -----------------------------
    > example_timeseries <- data.frame(
    +     time_s = t,
    +     random_vibration_g = random_vibration_scaled,
    +     sine_sweep = sine_sweep,
    +     inst_freq_hz = inst_freq
    + )
    >
    > random_psd_profile <- data.frame(
    +     freq_hz = f_bp,
    +     psd_g2_per_hz = psd_bp
    + )
    >
    > random_psd_estimate <- data.frame(
    +     freq_hz = f_est,
    +     psd_est_g2_per_hz = pxx
    + )
    >
    > # -----------------------------
    > # Save data to CSV
    > # -----------------------------
    > write.csv(example_timeseries,
    +           file = "d:/csv/example_timeseries.csv",
    +           row.names = FALSE)
    >
    > write.csv(random_psd_profile,
    +           file = "d:/csv/random_psd_profile.csv",
    +           row.names = FALSE)
    >
    > write.csv(random_psd_estimate,
    +           file = "d:/csv/random_psd_estimate.csv",
    +           row.names = FALSE)
    >
    > cat("\nData saved to CSV files in d:/csv/\n")
    >
    > # -----------------------------
    > # Create plots using ggplot2
    > # -----------------------------
    >
    > # Plot 1: Random Vibration PSD
    > plot1 <- ggplot() +
    +     geom_point(data = random_psd_profile,
    +                aes(x = freq_hz, y = psd_g2_per_hz, color = "Target PSD"),
    +                size = 3) +
    +     geom_line(data = random_psd_profile,
    +               aes(x = freq_hz, y = psd_g2_per_hz, color = "Target PSD"),
    +               linewidth = 1) +
    +     geom_line(data = random_psd_estimate,
    +               aes(x = freq_hz, y = psd_est_g2_per_hz, color = "Welch estimate"),
    +               alpha = 0.8, linewidth = 1) +
    +     scale_x_log10(breaks = c(1, 10, 100, 1000),
    +                   labels = c("1", "10", "100", "1000")) +
    +     scale_y_log10() +
    +     labs(x = "Frequency (Hz)",
    +          y = "PSD (gÂ²/Hz)",
    +          title = paste("Random Vibration PSD",
    +                        "target GRMS =", round(target_grms, 3),
    +                        "estimated GRMS =", round(estimated_grms, 3)),
    +          color = "") +
    +     theme_minimal() +
    +     theme(legend.position = c(0.85, 0.85),
    +           panel.grid.minor = element_line(linetype = "dotted")) +
    +     scale_color_manual(values = c("Target PSD" = "blue",
    +                                   "Welch estimate" = "red"))
    >
    > # Save plot
    > ggsave("d:/png/random_vibration_psd.png",
    +        plot = plot1, width = 8, height = 4.5, dpi = 160)
    Warning message:
    In scale_x_log10(breaks = c(1, 10, 100, 1000), labels = c("1", "10",  :
      log-10 transformation introduced infinite values.
    >
    > # Plot 2: Sine Sweep Signal
    > plot2 <- ggplot(example_timeseries, aes(x = time_s, y = sine_sweep)) +
    +     geom_line(color = "red", linewidth = 0.5) +
    +     labs(x = "Time (s)",
    +          y = "Amplitude",
    +          title = "Exponential Sine Sweep") +
    +     theme_minimal() +
    +     theme(panel.grid.minor = element_line(linetype = "dotted"))
    >
    > ggsave("d:/png/sine_sweep_signal.png",
    +        plot = plot2, width = 8, height = 4.5, dpi = 160)
    >
    > # Plot 3: Instantaneous Frequency
    > plot3 <- ggplot(example_timeseries, aes(x = time_s, y = inst_freq_hz)) +
    +     geom_line(color = "darkred", linewidth = 1) +
    +     labs(x = "Time (s)",
    +          y = "Frequency (Hz)",
    +          title = "Sine Sweep Instantaneous Frequency") +
    +     theme_minimal() +
    +     theme(panel.grid.minor = element_line(linetype = "dotted"))
    >
    > ggsave("d:/png/sine_sweep_frequency.png",
    +        plot = plot3, width = 8, height = 4.5, dpi = 160)
    >
    > # Plot 4: Random Vibration Time History (subset)
    > plot4 <- ggplot(example_timeseries[1:1000, ],
    +                 aes(x = time_s, y = random_vibration_g)) +
    +     geom_line(color = "blue", linewidth = 0.5) +
    +     labs(x = "Time (s)",
    +          y = "Acceleration (g)",
    +          title = "Random Vibration Time History (First 1000 samples)") +
    +     theme_minimal() +
    +     theme(panel.grid.minor = element_line(linetype = "dotted"))
    >
    > ggsave("d:/png/random_vibration_time_history.png",
    +        plot = plot4, width = 8, height = 4.5, dpi = 160)
    >
    > # -----------------------------
    > # Print summary statistics
    > # -----------------------------
    > cat("\n========================================\n")
    > cat("SUMMARY STATISTICS\n")
    > cat("========================================\n")
    >
    > cat("\nRandom Vibration Statistics:\n")
    > cat("Mean:", mean(random_vibration_scaled), "\n")
    > cat("Std Dev:", sd(random_vibration_scaled), "\n")
    > cat("Min:", min(random_vibration_scaled), "\n")
    > cat("Max:", max(random_vibration_scaled), "\n")
    > cat("GRMS (target):", target_grms, "\n")
    > cat("GRMS (estimated):", estimated_grms, "\n")
    >
    > cat("\nSine Sweep Statistics:\n")
    > cat("Mean:", mean(sine_sweep), "\n")
    > cat("Std Dev:", sd(sine_sweep), "\n")
    > cat("Min:", min(sine_sweep), "\n")
    > cat("Max:", max(sine_sweep), "\n")
    > cat("Frequency range:", min(inst_freq), "to", max(inst_freq), "Hz\n")
    >
    > cat("\n========================================\n")
    > cat("Analysis Complete!\n")
    > cat("========================================\n")

    NOTE: Processing of R statements complete

    319       run;
    NOTE: Procedure r step took :
          real time : 6.417
          cpu time  : 0.015



    NOTE: Submitted statements took :
          real time : 6.559
          cpu time  : 0.125

    /*              _
      ___ _ __   __| |
     / _ \ `_ \ / _` |
    |  __/ | | | (_| |
     \___|_| |_|\__,_|

    */
