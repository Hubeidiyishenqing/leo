const fs = require("fs");
const path = require("path");
const {
  Document, Packer, Paragraph, TextRun, Table, TableRow, TableCell,
  ImageRun, Header, Footer, AlignmentType, HeadingLevel,
  BorderStyle, WidthType, ShadingType, PageNumber, PageBreak,
  LevelFormat, TabStopType, TabStopPosition
} = require("docx");

// Read all 7 figures
const figDir = __dirname;
const fig1 = fs.readFileSync(path.join(figDir, "fig1_three_scheme_ber.png"));
const fig2 = fs.readFileSync(path.join(figDir, "fig2_ber_vs_elevation.png"));
const fig3 = fs.readFileSync(path.join(figDir, "fig3_ber_vs_scenario.png"));
const fig4 = fs.readFileSync(path.join(figDir, "fig4_overhead_vs_precomp.png"));
const fig5 = fs.readFileSync(path.join(figDir, "fig5_nav_accuracy_snr.png"));
const fig6 = fs.readFileSync(path.join(figDir, "fig6_nav_vs_elevation.png"));
const fig7 = fs.readFileSync(path.join(figDir, "fig7_dd_grid_comparison.png"));

// Helper functions
const border = { style: BorderStyle.SINGLE, size: 1, color: "CCCCCC" };
const borders = { top: border, bottom: border, left: border, right: border };
const headerBorder = { style: BorderStyle.SINGLE, size: 1, color: "2E75B6" };
const headerBorders = { top: headerBorder, bottom: headerBorder, left: headerBorder, right: headerBorder };

function makeCell(text, width, opts = {}) {
  const isHeader = opts.header || false;
  const bold = opts.bold || isHeader;
  const fontSize = opts.fontSize || 20;
  return new TableCell({
    borders: isHeader ? headerBorders : borders,
    width: { size: width, type: WidthType.DXA },
    shading: isHeader
      ? { fill: "2E75B6", type: ShadingType.CLEAR }
      : (opts.shading ? { fill: opts.shading, type: ShadingType.CLEAR } : undefined),
    margins: { top: 60, bottom: 60, left: 100, right: 100 },
    children: [
      new Paragraph({
        alignment: opts.align || AlignmentType.LEFT,
        children: [
          new TextRun({
            text: text,
            bold: bold,
            font: "Arial",
            size: fontSize,
            color: isHeader ? "FFFFFF" : (opts.color || "333333"),
          }),
        ],
      }),
    ],
  });
}

function makeParamRow(param, value, shading) {
  return new TableRow({
    children: [
      makeCell(param, 4000, { bold: true, shading }),
      makeCell(value, 5360, { shading }),
    ],
  });
}

function sectionTitle(text) {
  return new Paragraph({
    heading: HeadingLevel.HEADING_1,
    spacing: { before: 360, after: 200 },
    children: [new TextRun({ text, bold: true, font: "Arial", size: 32, color: "1A3C6E" })],
  });
}

function subTitle(text) {
  return new Paragraph({
    heading: HeadingLevel.HEADING_2,
    spacing: { before: 240, after: 120 },
    children: [new TextRun({ text, bold: true, font: "Arial", size: 26, color: "2E75B6" })],
  });
}

function bodyText(text) {
  return new Paragraph({
    spacing: { after: 120, line: 276 },
    children: [new TextRun({ text, font: "Arial", size: 21, color: "333333" })],
  });
}

function figureImage(imgData, w, h, caption) {
  return [
    new Paragraph({
      alignment: AlignmentType.CENTER,
      spacing: { before: 200, after: 80 },
      children: [
        new ImageRun({
          type: "png",
          data: imgData,
          transformation: { width: w, height: h },
          altText: { title: caption, description: caption, name: caption },
        }),
      ],
    }),
    new Paragraph({
      alignment: AlignmentType.CENTER,
      spacing: { after: 240 },
      children: [
        new TextRun({ text: caption, font: "Arial", size: 18, italics: true, color: "555555" }),
      ],
    }),
  ];
}

function bulletItem(text, ref) {
  return new Paragraph({
    numbering: { reference: ref, level: 0 },
    spacing: { after: 80, line: 276 },
    children: [new TextRun({ text, font: "Arial", size: 21, color: "333333" })],
  });
}

// Build document
const doc = new Document({
  styles: {
    default: { document: { run: { font: "Arial", size: 22 } } },
    paragraphStyles: [
      {
        id: "Heading1", name: "Heading 1", basedOn: "Normal", next: "Normal", quickFormat: true,
        run: { size: 32, bold: true, font: "Arial", color: "1A3C6E" },
        paragraph: { spacing: { before: 360, after: 200 }, outlineLevel: 0 },
      },
      {
        id: "Heading2", name: "Heading 2", basedOn: "Normal", next: "Normal", quickFormat: true,
        run: { size: 26, bold: true, font: "Arial", color: "2E75B6" },
        paragraph: { spacing: { before: 240, after: 120 }, outlineLevel: 1 },
      },
    ],
  },
  numbering: {
    config: [
      {
        reference: "bullets",
        levels: [{
          level: 0, format: LevelFormat.BULLET, text: "\u2022", alignment: AlignmentType.LEFT,
          style: { paragraph: { indent: { left: 720, hanging: 360 } } },
        }],
      },
      {
        reference: "numbers",
        levels: [{
          level: 0, format: LevelFormat.DECIMAL, text: "%1.", alignment: AlignmentType.LEFT,
          style: { paragraph: { indent: { left: 720, hanging: 360 } } },
        }],
      },
      {
        reference: "innovation",
        levels: [{
          level: 0, format: LevelFormat.DECIMAL, text: "%1.", alignment: AlignmentType.LEFT,
          style: { paragraph: { indent: { left: 720, hanging: 360 } } },
        }],
      },
    ],
  },
  sections: [
    {
      properties: {
        page: {
          size: { width: 12240, height: 15840 },
          margin: { top: 1440, right: 1200, bottom: 1440, left: 1200 },
        },
      },
      headers: {
        default: new Header({
          children: [
            new Paragraph({
              border: { bottom: { style: BorderStyle.SINGLE, size: 6, color: "2E75B6", space: 1 } },
              children: [
                new TextRun({ text: "LEO-NTN OTFS Simulation Progress Report", font: "Arial", size: 18, color: "2E75B6", bold: true }),
                new TextRun({ text: "\t2026-03-21", font: "Arial", size: 18, color: "888888" }),
              ],
              tabStops: [{ type: TabStopType.RIGHT, position: TabStopPosition.MAX }],
            }),
          ],
        }),
      },
      footers: {
        default: new Footer({
          children: [
            new Paragraph({
              alignment: AlignmentType.CENTER,
              children: [
                new TextRun({ text: "Page ", font: "Arial", size: 18, color: "888888" }),
                new TextRun({ children: [PageNumber.CURRENT], font: "Arial", size: 18, color: "888888" }),
              ],
            }),
          ],
        }),
      },
      children: [
        // ===== COVER =====
        new Paragraph({ spacing: { before: 2400 } }),
        new Paragraph({
          alignment: AlignmentType.CENTER,
          children: [new TextRun({ text: "LEO-NTN OTFS", font: "Arial", size: 52, bold: true, color: "1A3C6E" })],
        }),
        new Paragraph({
          alignment: AlignmentType.CENTER,
          spacing: { after: 120 },
          children: [new TextRun({ text: "Integrated Communication & Navigation", font: "Arial", size: 40, bold: true, color: "2E75B6" })],
        }),
        new Paragraph({
          alignment: AlignmentType.CENTER,
          spacing: { after: 600 },
          children: [new TextRun({ text: "Simulation Progress Report", font: "Arial", size: 36, color: "555555" })],
        }),
        new Paragraph({
          alignment: AlignmentType.CENTER,
          border: { top: { style: BorderStyle.SINGLE, size: 4, color: "2E75B6", space: 8 }, bottom: { style: BorderStyle.SINGLE, size: 4, color: "2E75B6", space: 8 } },
          spacing: { before: 200, after: 200 },
          children: [
            new TextRun({ text: "Version 5.0  |  March 21, 2026", font: "Arial", size: 24, color: "555555" }),
          ],
        }),
        new Paragraph({
          alignment: AlignmentType.CENTER,
          spacing: { before: 400 },
          children: [new TextRun({ text: "3GPP TR 38.811 NTN-TDL Channel  |  Shadowed-Rician LoS Fading", font: "Arial", size: 20, color: "888888" })],
        }),
        new Paragraph({
          alignment: AlignmentType.CENTER,
          children: [new TextRun({ text: "ICI-Aware DD-Domain Channel  |  2D Dirichlet Kernel Spreading", font: "Arial", size: 20, color: "888888" })],
        }),
        new Paragraph({
          alignment: AlignmentType.CENTER,
          children: [new TextRun({ text: "Quinn Fractional Estimator  |  LDPC (DVB-S.2, Rate 1/2)", font: "Arial", size: 20, color: "888888" })],
        }),

        // ===== PAGE BREAK =====
        new Paragraph({ children: [new PageBreak()] }),

        // ===== 1. RESEARCH OVERVIEW =====
        sectionTitle("1. Research Overview"),
        bodyText("This report presents the latest simulation results for the LEO-NTN OTFS integrated communication and navigation (ISAC) system. The simulation framework implements a complete OTFS transceiver chain with 3GPP-compliant NTN channel modeling, targeting journal publication with three key innovation points."),

        subTitle("1.1 Three Innovation Points"),
        new Paragraph({
          numbering: { reference: "innovation", level: 0 },
          spacing: { after: 80, line: 276 },
          children: [
            new TextRun({ text: "OTFS ISAC via DD-Domain Pilot: ", bold: true, font: "Arial", size: 21, color: "1A3C6E" }),
            new TextRun({ text: "Zero-overhead navigation by reusing the DD-domain pilot for both channel estimation and delay/Doppler sensing. Quinn fractional estimator enables sub-bin accuracy for both velocity and range estimation.", font: "Arial", size: 21, color: "333333" }),
          ],
        }),
        new Paragraph({
          numbering: { reference: "innovation", level: 0 },
          spacing: { after: 80, line: 276 },
          children: [
            new TextRun({ text: "DD-Domain Pilot + Doppler Pre-compensation Joint Design: ", bold: true, font: "Arial", size: 21, color: "1A3C6E" }),
            new TextRun({ text: "Satellite ephemeris-based bulk Doppler removal reduces the required guard band from 63 bins (infeasible, >100% overhead) to just a few bins (~0.7% overhead), making DD-domain pilot practical for LEO.", font: "Arial", size: 21, color: "333333" }),
          ],
        }),
        new Paragraph({
          numbering: { reference: "innovation", level: 0 },
          spacing: { after: 80, line: 276 },
          children: [
            new TextRun({ text: "OFDM/OTFS/OTFS-Pilot Three-Scheme Comparison: ", bold: true, font: "Arial", size: 21, color: "1A3C6E" }),
            new TextRun({ text: "Comprehensive performance comparison under realistic LEO NTN conditions with ICI-aware channel modeling, demonstrating OTFS-Pilot superiority in both communication BER and navigation accuracy.", font: "Arial", size: 21, color: "333333" }),
          ],
        }),

        // ===== 2. SYSTEM PARAMETERS =====
        new Paragraph({ children: [new PageBreak()] }),
        sectionTitle("2. System Parameters"),
        bodyText("The simulation uses the following system configuration, aligned with 3GPP NR NTN specifications:"),

        new Table({
          width: { size: 9360, type: WidthType.DXA },
          columnWidths: [4000, 5360],
          rows: [
            new TableRow({
              children: [
                makeCell("Parameter", 4000, { header: true }),
                makeCell("Value", 5360, { header: true }),
              ],
            }),
            makeParamRow("Carrier Frequency (fc)", "2 GHz (S-band)", "F0F5FA"),
            makeParamRow("System Bandwidth", "100 MHz"),
            makeParamRow("Subcarrier Spacing (SCS)", "120 kHz", "F0F5FA"),
            makeParamRow("Delay Bins (N)", "1024"),
            makeParamRow("Doppler Bins (M)", "128", "F0F5FA"),
            makeParamRow("CP Ratio", "7%"),
            makeParamRow("Modulation", "16-QAM", "F0F5FA"),
            makeParamRow("Channel Coding", "LDPC DVB-S.2, Rate 1/2"),
            makeParamRow("UE Velocity", "26,000 km/h (LEO orbital)", "F0F5FA"),
            makeParamRow("Max Doppler Shift (fd)", "~48,148 Hz"),
            makeParamRow("Delay Resolution", "~8.14 ns/bin", "F0F5FA"),
            makeParamRow("Doppler Resolution", "~0.877 Hz/bin"),
            makeParamRow("Satellite Altitude", "600 km (LEO)", "F0F5FA"),
            makeParamRow("Channel Model", "3GPP TR 38.811 NTN-TDL-C"),
            makeParamRow("LoS Fading", "Shadowed-Rician (Loo model)", "F0F5FA"),
            makeParamRow("Shadow Std (ITU-R P.681)", "Elevation & scenario dependent"),
            makeParamRow("DD Channel Application", "2D Dirichlet kernel spreading", "F0F5FA"),
            makeParamRow("Sensing Estimator", "Quinn fractional (delay + Doppler)"),
            makeParamRow("OFDM Equalizer", "TF-MMSE with post-AFC channel", "F0F5FA"),
            makeParamRow("OTFS Equalizer", "TF-MMSE (conventional) / DD-MP (pilot)"),
          ],
        }),

        // ===== 3. SIMULATION RESULTS =====
        new Paragraph({ children: [new PageBreak()] }),
        sectionTitle("3. Simulation Results"),
        bodyText("This section presents the 7 simulation figures covering communication BER performance, pilot overhead analysis, navigation accuracy, and DD-domain visualization."),

        // --- Fig 1 ---
        subTitle("3.1 Three-Scheme BER vs Eb/N0 (Fig. 1)"),
        ...figureImage(fig1, 550, 420, "Fig. 1: Three-Scheme BER vs Eb/N0 (NTN-TDL-C, Suburban, Elev=50\u00B0)"),
        bodyText("Figure 1 compares the BER performance of three schemes across Eb/N0. The OTFS-Pilot scheme with DD-domain Message Passing (MP) equalization achieves the best performance, reaching BER < 10\u207B\u2074 at moderate SNR. The coded OTFS-Pilot (C-OTFS-P) with LDPC further extends the coding gain. OFDM and conventional OTFS using TF-MMSE equalization exhibit high error floors due to the severe ICI from LEO Doppler spread (~48 kHz) that cannot be fully equalized in the TF domain. This demonstrates the fundamental advantage of DD-domain processing for high-mobility LEO channels."),

        // --- Fig 2 ---
        new Paragraph({ children: [new PageBreak()] }),
        subTitle("3.2 BER vs Satellite Elevation Angle (Fig. 2)"),
        ...figureImage(fig2, 600, 300, "Fig. 2: BER vs Elevation Angle (NTN-TDL-C, Suburban)"),
        bodyText("Figure 2 shows BER variation with satellite elevation angle at two SNR operating points (12 dB and 16 dB). The OTFS-Pilot scheme shows improving BER with higher elevation angles, consistent with the elevation-dependent K-factor increase in the NTN-TDL-C profile. At 16 dB, OTFS-Pilot BER improves from ~5\u00D710\u207B\u00B3 at 10\u00B0 to ~3\u00D710\u207B\u2074 at 90\u00B0. OFDM and conventional OTFS remain at high BER (~0.4\u20130.5) regardless of elevation, confirming the ICI-limited error floor is independent of channel K-factor."),

        // --- Fig 3 ---
        subTitle("3.3 BER vs Propagation Scenario (Fig. 3)"),
        ...figureImage(fig3, 580, 380, "Fig. 3: BER by Scenario (Eb/N0=15 dB, Elev=50\u00B0)"),
        bodyText("Figure 3 compares BER across four propagation scenarios at fixed Eb/N0=15 dB and 50\u00B0 elevation. OTFS-Pilot (DD-MP) achieves BER ranging from 1.9\u00D710\u207B\u00B2 (Dense Urban) to 5.4\u00D710\u207B\u00B3 (Rural), reflecting the scenario-dependent Shadowed-Rician fading and K-factor variation. Rural environments with higher K-factor and lower shadow fading std yield the best performance. OFDM and OTFS (TF-MMSE) remain at ~0.45\u20130.50 across all scenarios."),

        // --- Fig 4 ---
        new Paragraph({ children: [new PageBreak()] }),
        sectionTitle("4. Pilot Design & Overhead Analysis"),
        subTitle("4.1 Pilot Overhead vs Pre-compensation Ratio (Fig. 4)"),
        ...figureImage(fig4, 560, 280, "Fig. 4: Pilot + Guard Overhead vs Doppler Pre-compensation Ratio"),
        bodyText("Figure 4 demonstrates the critical role of Doppler pre-compensation in enabling DD-domain pilot placement. Without pre-compensation, the required Doppler guard band is k_Guard=63 bins (exceeding M/2=64), making the pilot pattern infeasible with >100% overhead. With increasing pre-compensation ratio, the guard band shrinks dramatically. At 95%+ pre-compensation (achievable with satellite ephemeris), the total pilot+guard overhead drops to ~0.7%, which is negligible compared to conventional pilot schemes."),
        bodyText("This result validates Innovation Point 2: Doppler pre-compensation is not merely an optimization but an absolute necessity for DD-domain pilot in LEO satellite channels."),

        // --- Fig 5 ---
        new Paragraph({ children: [new PageBreak()] }),
        sectionTitle("5. Navigation Accuracy Results"),
        subTitle("5.1 Velocity & Range RMSE vs SNR (Fig. 5)"),
        ...figureImage(fig5, 560, 280, "Fig. 5: Navigation Accuracy (Velocity & Range RMSE) vs Eb/N0"),
        bodyText("Figure 5 evaluates the navigation sensing performance of OTFS DD-domain pilot versus OFDM pilot (SFFT+Quinn) as a function of SNR. The OTFS DD pilot with fractional 2D Dirichlet kernel spreading achieves superior velocity and range estimation accuracy, with RMSE decreasing steadily with SNR. The CRB (Cramer-Rao Bound) is shown as the theoretical lower bound. The OFDM SFFT+Quinn approach shows less stable performance, with higher variance across SNR points due to the ICI-corrupted TF-domain observations degrading the 2D-FFT peak quality."),

        // --- Fig 6 ---
        subTitle("5.2 Navigation Accuracy vs Elevation Angle (Fig. 6)"),
        ...figureImage(fig6, 520, 420, "Fig. 6: Navigation RMSE & OTFS Advantage Factor vs Elevation Angle"),
        bodyText("Figure 6 presents a four-panel analysis of navigation accuracy versus elevation angle at Eb/N0=15 dB. Panels (a) and (b) show velocity and range RMSE for both schemes. Panels (c) and (d) show the OTFS advantage factor (ratio of OFDM RMSE to OTFS RMSE). The OTFS DD pilot demonstrates up to 7\u00D7 velocity advantage and up to 5.1\u00D7 range advantage at certain elevation angles. The advantage is most pronounced at higher elevations where the stronger LoS component provides better DD-domain channel estimation conditions."),

        // --- Fig 7 ---
        new Paragraph({ children: [new PageBreak()] }),
        sectionTitle("6. DD Grid Visualization"),
        subTitle("6.1 DD Grid Pilot Pattern (Fig. 7)"),
        ...figureImage(fig7, 600, 300, "Fig. 7: DD Grid - No Pre-comp (Infeasible) vs With Pre-comp (0.7% Overhead)"),
        bodyText("Figure 7 provides a visual comparison of the DD-domain pilot placement with and without Doppler pre-compensation. Panel (a) shows the case without pre-compensation: the required Doppler guard band k_Guard=63 exceeds M/2=64, making the entire grid infeasible (shown in orange with INFEASIBLE warning). Panel (b) shows the result with pre-compensation: the pilot and guard region (highlighted with cyan dashed border) occupies only a small fraction of the grid (0.7% overhead), demonstrating the practical feasibility of DD-domain ISAC in LEO channels."),

        // ===== 7. KEY TECHNICAL FEATURES =====
        new Paragraph({ children: [new PageBreak()] }),
        sectionTitle("7. Key Technical Features"),

        subTitle("7.1 ICI-Aware Channel Modeling"),
        bodyText("The simulation applies the multipath channel directly in the Delay-Doppler domain via 2D circular convolution (applyChannelDD), capturing the full Inter-Carrier Interference (ICI) caused by high Doppler spread. This is critical for LEO scenarios where fd/SCS \u2248 40%, far exceeding the quasi-static assumption of conventional TF-domain channel multiplication. The TF-MMSE equalizer uses a post-AFC channel matrix (buildTF_AFC) that accounts for residual Doppler after bulk pre-compensation."),

        subTitle("7.2 Shadowed-Rician LoS Fading"),
        bodyText("The LoS component uses the Shadowed-Rician (Loo) model where the direct-path amplitude follows lognormal fading, capturing ionospheric scintillation and tropospheric shadowing effects on satellite channels. Shadow fading standard deviation is elevation-angle and scenario-dependent, following ITU-R P.681 tables. This provides realistic trial-to-trial variation absent from deterministic Rician models."),

        subTitle("7.3 2D Dirichlet Kernel Spreading"),
        bodyText("For navigation accuracy evaluation, the DD-domain channel application uses fractional mode with 2D Dirichlet kernel spreading coefficients for both delay and Doppler dimensions. This models Inter-Delay Interference (IDI) and Inter-Doppler Interference enabling the Quinn fractional estimator to achieve sub-bin accuracy. The spreading half-width is Ni=2 in both dimensions, resulting in a (2\u00D72+1)\u00B2 = 25 coefficient grid per path."),

        subTitle("7.4 Quinn Fractional Estimator"),
        bodyText("Both OFDM (SFFT+Quinn) and OTFS (DD pilot) sensing employ the Quinn fractional correction algorithm for sub-bin delay and Doppler estimation. This ensures a fair comparison between schemes. The estimator extracts the fractional offset from the ratio of adjacent spectral bins around the detected peak, providing interpolation accuracy significantly beyond the grid resolution."),

        // ===== 8. CURRENT STATUS =====
        sectionTitle("8. Current Status & Observations"),
        bulletItem("OTFS-Pilot (DD-MP) demonstrates strong BER performance with proper error floor behavior across all scenarios and elevation angles.", "bullets"),
        bulletItem("OFDM and conventional OTFS with TF-MMSE show ICI-limited error floors (~0.45), confirming that TF-domain equalization is insufficient for extreme LEO Doppler.", "bullets"),
        bulletItem("Doppler pre-compensation reduces pilot overhead from infeasible (>100%) to negligible (~0.7%).", "bullets"),
        bulletItem("OTFS DD pilot achieves up to 7\u00D7 velocity and 5.1\u00D7 range advantage over OFDM SFFT+Quinn sensing.", "bullets"),
        bulletItem("Shadowed-Rician fading introduces realistic LoS variation, with Rural scenarios showing best performance due to low shadow std.", "bullets"),
        bulletItem("Navigation accuracy improves with SNR and shows elevation-dependent behavior consistent with K-factor variation.", "bullets"),

        // ===== 9. FIGURE INDEX =====
        new Paragraph({ children: [new PageBreak()] }),
        sectionTitle("9. Figure Index"),
        new Table({
          width: { size: 9840, type: WidthType.DXA },
          columnWidths: [1200, 5640, 3000],
          rows: [
            new TableRow({
              children: [
                makeCell("Figure", 1200, { header: true, align: AlignmentType.CENTER }),
                makeCell("Description", 5640, { header: true }),
                makeCell("Innovation", 3000, { header: true, align: AlignmentType.CENTER }),
              ],
            }),
            new TableRow({ children: [
              makeCell("Fig. 1", 1200, { align: AlignmentType.CENTER, shading: "F0F5FA" }),
              makeCell("Three-Scheme BER vs Eb/N0", 5640, { shading: "F0F5FA" }),
              makeCell("Innovation 3", 3000, { align: AlignmentType.CENTER, shading: "F0F5FA" }),
            ]}),
            new TableRow({ children: [
              makeCell("Fig. 2", 1200, { align: AlignmentType.CENTER }),
              makeCell("BER vs Satellite Elevation Angle", 5640),
              makeCell("Innovation 3", 3000, { align: AlignmentType.CENTER }),
            ]}),
            new TableRow({ children: [
              makeCell("Fig. 3", 1200, { align: AlignmentType.CENTER, shading: "F0F5FA" }),
              makeCell("BER vs Propagation Scenario", 5640, { shading: "F0F5FA" }),
              makeCell("Innovation 3", 3000, { align: AlignmentType.CENTER, shading: "F0F5FA" }),
            ]}),
            new TableRow({ children: [
              makeCell("Fig. 4", 1200, { align: AlignmentType.CENTER }),
              makeCell("Pilot Overhead vs Pre-compensation Ratio", 5640),
              makeCell("Innovation 2", 3000, { align: AlignmentType.CENTER }),
            ]}),
            new TableRow({ children: [
              makeCell("Fig. 5", 1200, { align: AlignmentType.CENTER, shading: "F0F5FA" }),
              makeCell("Navigation Accuracy (Velocity & Range RMSE) vs SNR", 5640, { shading: "F0F5FA" }),
              makeCell("Innovation 1", 3000, { align: AlignmentType.CENTER, shading: "F0F5FA" }),
            ]}),
            new TableRow({ children: [
              makeCell("Fig. 6", 1200, { align: AlignmentType.CENTER }),
              makeCell("Navigation Accuracy vs Elevation Angle + Advantage Factor", 5640),
              makeCell("Innovation 1", 3000, { align: AlignmentType.CENTER }),
            ]}),
            new TableRow({ children: [
              makeCell("Fig. 7", 1200, { align: AlignmentType.CENTER, shading: "F0F5FA" }),
              makeCell("DD Grid Pilot Pattern Visualization", 5640, { shading: "F0F5FA" }),
              makeCell("Innovation 2", 3000, { align: AlignmentType.CENTER, shading: "F0F5FA" }),
            ]}),
          ],
        }),
      ],
    },
  ],
});

// Generate
const outPath = path.join(figDir, "LEO_OTFS_Report_v5.docx");
Packer.toBuffer(doc).then((buffer) => {
  fs.writeFileSync(outPath, buffer);
  console.log("Report generated: " + outPath);
});
