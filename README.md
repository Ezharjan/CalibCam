# CalibCam — manual camera calibration in the browser

CalibCam is a single-page web tool for calibrating a camera **by hand** from one image: a public webcam snapshot URL, an image or video file, or your own webcam. You annotate geometry you can see — straight edges, parallel lines, points with known world coordinates — and the tool estimates the camera with classical, verifiable geometry. **No machine learning / deep learning is used anywhere.**

Output is a standard OpenCV-style calibration:

| Quantity | Meaning |
|---|---|
| `K` = `[[fx,0,cx],[0,fy,cy],[0,0,1]]` | intrinsics: focal length (px) and principal point (px) |
| `dist` = `[k1, k2, p1, p2, k3]` | OpenCV radial–tangential lens distortion |
| `R`, `rvec`, `t` | extrinsics, world → camera: `x_cam = R·X_world + t` |
| `C = −Rᵀt` | camera centre in world coordinates |
| `H` | ground-plane homography, plane `(X,Y,1)` → undistorted pixel |

---

## Contents

- [CalibCam — manual camera calibration in the browser](#calibcam--manual-camera-calibration-in-the-browser)
  - [Contents](#contents)
  - [1. Files \& deployment](#1-files--deployment)
  - [2. Conventions (read this once)](#2-conventions-read-this-once)
  - [3. The interface](#3-the-interface)
  - [4. Loading a frame](#4-loading-a-frame)
    - [4.1 A webcam / image URL (main input)](#41-a-webcam--image-url-main-input)
    - [4.2 Many URLs (batch annotation)](#42-many-urls-batch-annotation)
    - [4.3 Other sources](#43-other-sources)
  - [5. Tutorial — calibrating a camera by hand](#5-tutorial--calibrating-a-camera-by-hand)
    - [Step 1 · Lens distortion (plumb lines)](#step-1--lens-distortion-plumb-lines)
    - [Step 2 · Vanishing points](#step-2--vanishing-points)
    - [Step 3 · Point correspondences](#step-3--point-correspondences)
    - [Step 4 · Ground-plane homography (optional)](#step-4--ground-plane-homography-optional)
    - [Step 5 · Verify](#step-5--verify)
    - [Step 6 · Export](#step-6--export)
  - [6. Recipes for typical scenes](#6-recipes-for-typical-scenes)
  - [7. Accuracy checklist](#7-accuracy-checklist)
  - [8. Reading the numbers](#8-reading-the-numbers)
  - [9. Practise with the synthetic test scene](#9-practise-with-the-synthetic-test-scene)
  - [10. Keyboard \& mouse reference](#10-keyboard--mouse-reference)
  - [11. Export format](#11-export-format)
  - [12. Using the result in OpenCV](#12-using-the-result-in-opencv)
  - [13. Methods (what the math does)](#13-methods-what-the-math-does)
  - [14. Troubleshooting](#14-troubleshooting)
  - [15. Privacy](#15-privacy)
  - [16. Limitations](#16-limitations)

---

## 1. Files & deployment

| File | Needed | Purpose |
|---|---|---|
| `Camera Calibration.dc.html` | **yes** | the app (rename to `index.html` for GitHub Pages) |
| `support.js` | **yes** | runtime the page needs |
| `calib.js` | **yes** | all calibration math (no dependencies) |
| `webcams.csv` | optional | a URL list that is auto-loaded if present — omit it to keep a list private (you can still import it at runtime) |

Everything else in an export (`screenshots/`, `uploads/`, `.thumbnail`) is not needed.

**GitHub Pages:** put the three required files in the repository root, rename the HTML file to `index.html`, enable Pages. Pages is served over `https://`, so your own webcam works there.

**Locally:** do **not** double-click the HTML file (browsers restrict `file://` pages). Serve the folder instead:

```bash
cd path/to/folder
python3 -m http.server 8000
# open http://localhost:8000/Camera%20Calibration.dc.html
```

`localhost` and `https://` count as secure contexts, which browsers require for camera access.

---

## 2. Conventions (read this once)

- **Pixel coordinates** follow OpenCV: the *centre* of the top-left pixel is `(0, 0)`, `u` grows right, `v` grows down. All coordinates you see are sub-pixel.
- **Camera frame**: x right, y down, z forward (OpenCV).
- **World frame**: you choose it. It must be **right-handed**. For a ground plane, a natural choice is X and Y on the ground and **Z pointing up**. If the tool reports a *reflection* or *left-handed* world, you have swapped an axis.
- **Units** of `t` and `C` are whatever units you type for world coordinates (m, mm, tiles, …). Pixel quantities are always in pixels.
- **Distortion coefficients are relative to `fx`.** Whenever the tool changes `fx` while holding distortion fixed, it rescales `k1,k2,k3,p1,p2` so that the *pixel-space* distortion is unchanged (`k1·s², k2·s⁴, k3·s⁶, p·s` with `s = f_new/f_old`). You never have to redo Lens because f changed.
- **The camera on the right is shared.** Every step reads it and, when you press its solve button, writes its result back. The *Camera* card always shows the current state and which step last changed it. **⌘/Ctrl-Z** undoes annotations *and* camera changes.

---

## 3. The interface

Two layouts, switch at the top right (the choice is remembered):

- **A · Workbench** — tool rail on the left, image in the middle, inspector on the right. For people who know the workflow and jump between tools.
- **B · Guided** — numbered stepper across the top with ✓ marks for completed steps, and Back / Next at the bottom of the inspector. Use this the first few times.

Common parts:

- **Canvas toolbar** (top-left of the image)
  - **View: raw / undistorted** — show the original image or the image rectified with the current lens model. You can annotate in either view; coordinates are always stored in raw pixels.
  - **Grid overlay** — projects a world grid (Z = 0), coloured axes (X red, Y green, Z blue) and an optional unit cube through the *full* camera model, distortion included.
  - **Residuals ×10** — for each point with world coordinates, draws a magenta vector from the clicked position towards its reprojection, **magnified 10×**, plus a small circle at the true reprojection.
  - **Loupe** — 12× magnifier that follows the cursor, with crosshair and sub-pixel readout. Turn it off if it covers what you need.
  - **Fit** — fit the image to the window (also **F**).
- **Status bar** (bottom) — cursor position in raw pixels with two decimals; on the Plane tool also the world plane coordinates under the cursor; zoom; a hint for the current tool.
- **Camera card** (bottom of the inspector) — `K`, `dist`, `rvec`, `t`, `C`, field of view; *Download JSON*, *Copy JSON*, *Reset*.

---

## 4. Loading a frame

### 4.1 A webcam / image URL (main input)

1. Open **Source**.
2. Paste the snapshot URL into **Webcam / image URL** (e.g. `https://…/latest-frame.jpg`) and press **Load** or Enter. The same field is on the empty canvas.
3. The frame appears; its name, size and time are shown in the Source panel. Recently used URLs are listed under the field for one-click reuse.
4. **Fetch latest frame (keep annotations)** (bottom of the canvas) re-downloads the same URL. Annotations stay as long as the image size doesn't change — useful for fixed cameras, or to wait for better light.

The URL must return a **still image** (JPEG/PNG/WebP). HTML pages, MJPEG streams and HLS/RTSP video are not supported — find the camera's snapshot endpoint (often `…/snapshot.jpg`, `…/latest.jpg`, `…/image.jpg`).

**About the image proxy.** Browsers only let a page read an image's pixels if the camera server allows it (CORS). Most public webcam servers don't — in testing, *every* host in a 254-camera list refused direct pixel access. The tool therefore tries the URL directly first and, if that fails, fetches it through an image proxy:

- **Fallback via image proxy** (toggle) — on by default.
- **Proxy template** — `{url}` is replaced with the URL-encoded camera URL. Leave it blank to use the public default `https://images.weserv.nl/?url={url}`, or enter your own proxy (see [Privacy](#15-privacy)).
- With the proxy **off**, only CORS-enabled `https://` cameras will load. An `http://` camera can never load directly from an `https://` page (mixed content).

The proxy does not change image geometry (no resizing or cropping), so calibration stays valid. It may re-encode the JPEG, which adds compression noise comparable to the camera's own.

### 4.2 Many URLs (batch annotation)

If you annotate many cameras, click **Annotating many cameras? Import a URL list (CSV)…**. Accepted formats:

```csv
ID,Name,ImageURL,Latitude,Longitude
lake_cam,Lake Cam,https://example.com/lake.jpg,20.30,-103.27
```

or simply one URL per line. A searchable **Imported URL list** appears; click a row to load it. Loading a *different* camera clears the annotations; refreshing the *same* camera keeps them. The list stays in your browser only. If you deploy a `webcams.csv` next to the HTML file, it is loaded automatically.

### 4.3 Other sources

- **Open image / video** or drag & drop. For video, scrub with the slider or step ±1/30 s; each frame of the same size keeps the annotations.
- **Webcam** — your own camera. **Capture frame** freezes a frame for annotation (sub-pixel snapping only works on captured frames); **Resume live** goes back to the feed. This requires a secure context and your permission (see [Troubleshooting](#14-troubleshooting)).
- **Synthetic test scene** — a ray-traced room corner with known ground truth, for practice (see [§9](#9-practise-with-the-synthetic-test-scene)).
- **Import calibration JSON** — reload a previous export (camera and annotations).

---

## 5. Tutorial — calibrating a camera by hand

The recommended order is **Lens → Vanishing points → Points → Verify → Export**. Each step improves on the previous one; skip the ones your scene can't support (see [§6](#6-recipes-for-typical-scenes)).

Before you start: zoom in (**mouse wheel**) for every click. A click at 27% zoom is worth far less than one at 300%. Pan with **right-drag** or **Space + drag**.

### Step 1 · Lens distortion (plumb lines)

Goal: estimate `k1` (and `k2`, `k3`) from edges that are **straight in the real world** but bowed in the image.

1. Open **Lens**. Keep **Snap to edges (sub-pixel)** on.
2. Find long straight structures: building edges, roof lines, poles, road markings, the horizon line of a flat sea, window frames, door frames.
3. Click points along one edge, **5–15 points**, spread end to end. Press **Enter** (or double-click) to finish the line. When snapping is on, each point is pulled perpendicular to the edge onto the strongest gradient, with parabolic sub-pixel interpolation.
4. Repeat for **3–8 lines**. What matters most:
   - **Lines near the image borders and corners.** Distortion grows with radius; lines through the centre carry almost no information.
   - **Different orientations** (horizontal, vertical, diagonal).
   - Avoid edges that are actually curved (cables, arches, wavy roofs) and the far horizon over hilly land.
5. Choose the **radial model**:
   - `k1` — most webcams and normal lenses; most robust.
   - `k1 k2` — wide-angle lenses, or when `k1` leaves visible curvature.
   - `k1 k2 k3` — fisheye-like lenses with many long border lines only.
6. Press **Solve distortion**. Read the result:
   - **Straightness RMS `before → after`** — mean perpendicular deviation of your points from a straight line, in pixels. A good result is **≤ 0.3–0.5 px** after.
   - Each coefficient with **± σ** (1-sigma from the LM covariance). If σ is as large as the value itself, the coefficient is not determined — use a simpler model or add border lines.
   - The **per-line deviation** in the list. A line far above the others is probably not straight in reality — delete it (×) and re-solve.
7. Check with **View: undistorted**: your traced edges should now look straight.

Notes: the distortion centre is the current principal point. **Zero k** resets distortion. Tangential distortion (`p1, p2`) is not estimated here (plumb lines constrain it poorly); it can be refined in Points.

### Step 2 · Vanishing points

Goal: get `f`, the principal point and the camera rotation `R` from sets of parallel world lines.

1. Open **Vanishing points**. Pick an **axis** (X / Y / Z buttons; colours red / green / blue).
2. Draw segments along edges **parallel to that world axis**: **drag** from one end to the other, or click both endpoints. With **Fit segments to edges** on, each segment is fitted to the image edge with sub-pixel accuracy — in undistorted space, so distortion from Step 1 is respected. Segments drawn after Lens appear slightly curved: that is the straight 3-D edge seen through the lens.
3. Draw **≥ 2 lines per axis, ideally 3–6**, for **two or three mutually orthogonal axes**. Typical choice: X and Y along the two horizontal directions of a building or room, Z vertical.
   - Use **long** segments.
   - Use lines that **converge strongly** (far apart in the image, receding in depth). Nearly-parallel image lines put the vanishing point near infinity and make `f` uncertain.
   - Spread lines across the image.
4. Choose the **principal point** mode:
   - **Image centre** — default; correct to within a few pixels for most cameras. Needs 2 axes.
   - **Keep current** — use `cx, cy` from an earlier step.
   - **From 3 VPs** — estimate it as the orthocentre of the three vanishing points. Needs 3 axes with strong perspective on all three; otherwise prefer Image centre.
5. Press **Solve f, c, R**. Read the result:
   - Each VP position (or ∞) and its **RMS endpoint distance** (px). With good lines this is **< 0.5 px**.
   - Per-line **angle residual** in the list (°). Remove lines with outlying values.
   - **f ± σ**: σ comes from 120 Monte-Carlo trials in which every endpoint is perturbed by the noise level measured from your own residuals. Aim for **σ/f < 1–2 %**.
   - Warnings: zero distortion (do Lens first), > 5 % uncertainty, axes far from orthogonal.
6. The world axes have a ±90° ambiguity around Z. If the overlay grid is rotated or mirrored with respect to what you intended, press **Rotate XY 90°** (up to 3 times).
7. Translation **cannot** be observed from vanishing points. The camera centre is kept (or placed 10 units in front of the origin the first time). Step 3 or the Verify sliders set it.

### Step 3 · Point correspondences

Goal: the most accurate step — refine everything jointly by minimising **reprojection error** over points whose world coordinates you know.

1. Open **Points**. Keep **Snap to corners (sub-pixel)** on: each click is refined to the nearest corner with the gradient-orthogonality method (same principle as OpenCV `cornerSubPix`). If there is no corner within ~5 px, the click is kept as is.
2. Click a feature, then type its **X Y Z** in the table (blank Z = 0). Good features: corners of road markings, tiles, windows, building corners, survey marks, posts at known positions. Get world coordinates from site plans, measurements, maps or known standard dimensions.
3. **Fine-tune**: select a point (click it or its row) and nudge with the arrow keys (**0.1 px**, Shift = 1 px), or drag it. Use the loupe.
4. How many points?
   - **≥ 6 non-coplanar** (some with Z ≠ 0) → full DLT start, can estimate everything.
   - **≥ 4 coplanar** (all on one plane) → homography start; can estimate `f` and pose but **not** the principal point (held fixed automatically).
   - More is better: **10–30** well-spread points. Cover the whole image, including the borders, if you refine distortion.
   - The line under the options shows `n pts · coplanar/3D · residuals / params` — keep residuals well above parameters (≥ 2–3×).
5. Choose what to estimate:
   - **Focal**: *Est. fx=fy* (square pixels, recommended), *Est. fx, fy*, or *Fixed*.
   - **Principal pt**: *Fixed* (recommended unless you have many 3-D points) or *Estimate*.
   - **Distortion**: *keep* (hold the Lens result, recommended), or refine `k1` / `k1k2` / `+p1p2` / *all 5* jointly. Only refine distortion with many points spread to the borders.
   - **Start from**: *Linear solution* (DLT/homography, recommended) or *Current camera* (use the VP or slider estimate as the start — needed when you have fewer than 6 non-coplanar points).
6. Toggle a point in or out of the solve by clicking its **number**. Delete with ×.
7. Press **Calibrate**. Read the result:
   - **Reprojection RMS (Euclidean / pt)** — `sqrt(mean(du² + dv²))`, the same definition as OpenCV. With sub-pixel corners and correct coordinates, **< 0.5 px** is good and **< 1 px** acceptable.
   - **Max error** and the per-point **err** column (green < 0.5 px, amber > 2 px). An amber point usually has a typo in its world coordinates or was clicked on the wrong feature. Fix it or disable it and re-solve.
   - **σ** of every estimated parameter, rotation (°) and `|t|`. Large σ means the data doesn't constrain that parameter — fix it or add points.
   - The log (initialisation method, iterations, degrees of freedom) and warnings (near fronto-parallel plane, low redundancy, points behind the camera).

### Step 4 · Ground-plane homography (optional)

Use this when you need a pixel ↔ ground mapping (traffic, sports, surveillance).

1. Open **Plane**. Only points with **Z = 0** are used; others are excluded automatically and counted.
2. Press **Solve homography** (≥ 4 points, no 3 collinear; 8+ recommended).
3. Check the **RMS error** and the purple plane grid drawn through `H`.
4. **Hover** over the image: the status bar shows the plane coordinates `(X, Y)` under the cursor, computed through `H⁻¹` from the undistorted pixel. Use this to spot-check known distances.

`H` maps plane `(X, Y, 1)` to **undistorted** pixels. To use it on raw pixels, undistort them first with `K` and `dist` (see [§12](#12-using-the-result-in-opencv)).

### Step 5 · Verify

1. Open **Verify**. The world grid (Z = 0, from the origin, `grid −/+` to change its size), axes and the optional **unit cube** are projected through the full model.
2. Check visually: grid lines should follow ground features, vertical cube edges should line up with real verticals, the X/Y/Z axes should point where you defined them. Check near the **image borders** too — that is where a wrong distortion model shows.
3. **Sliders**: `f`, `cx`, `cy`, `k1`, `k2`, the rotation angles α β γ (`R = Rz(γ)·Ry(β)·Rx(α)`, world→camera) and the camera centre `C`. Drag for coarse changes; type exact values in the fields (Enter applies). The live **points RMS** updates as you move. Use the sliders to
   - seed a pose when you have too few points for a linear solution, then return to Points with *Start from: Current camera*;
   - do a fully manual calibration when the scene has no measurable points.
4. With the synthetic scene, a **ground-truth table** shows the error of every parameter.

### Step 6 · Export

**Export JSON** (header), **Download JSON** or **Copy JSON** (Camera card). The file contains the camera, the per-step results, all annotations and a description of the conventions, so it can be re-imported with **Import calibration JSON** and the work continued.

---

## 6. Recipes for typical scenes

**Street / building webcam (most common)**
Lens on building edges and poles near the borders → Vanishing points with X/Y along the two façade directions and Z on verticals (Image centre) → Points on road-marking corners and building corners measured from a site plan or map (focal *Est. fx=fy*, pp *Fixed*, distortion *keep*) → Verify.

**Flat ground only (sports field, car park, plaza)**
Lens if any straight long edges exist → Points with ≥ 8 coplanar points on known markings (Z = 0) → `f` and pose from the plane (keep the principal point fixed, avoid a near fronto-parallel view) → Plane for the homography.

**Landscape / mountain camera (no man-made geometry)**
There is usually not enough geometry for full calibration. Use the horizon only if it is a true sea horizon. Set `f` from the camera's known field of view if available (Verify sliders), then Points with known landmarks (peaks, towers) using coordinates converted to a local metric frame, *Start from: Current camera*. Expect larger uncertainties and read the σ values.

**Indoor room / corridor**
Ideal for three vanishing points (two wall directions and the vertical) → *From 3 VPs* principal point is feasible → Points on tiles or door corners.

**Your own webcam / phone**
Print or display a checkerboard with known square size, capture a frame with the board **tilted ~30–50°** and covering much of the image → Lens on the board edges → Points on board corners (Z = 0) → Calibrate.

---

## 7. Accuracy checklist

- [ ] Every click made at **≥ 200 % zoom**, snapping on, loupe checked.
- [ ] Lens solved **before** Vanishing points (otherwise `f` can be biased by 5–10 %).
- [ ] Plumb lines cover the **borders and corners**, not the centre.
- [ ] VP lines are **long**, **spread out** and **strongly converging**; ≥ 3 per axis.
- [ ] World coordinates are in **one right-handed frame** and **one unit**, double-checked for typos.
- [ ] Points **spread over the whole image**, with depth variation (or a tilted plane).
- [ ] Residuals ≥ 2–3 × parameters; no parameter estimated that the data can't support.
- [ ] Reprojection RMS **< 0.5–1 px**, no outlier point, σ values small.
- [ ] Overlay checked **near the borders** and with the unit cube.
- [ ] For a fixed camera: **Fetch latest frame** at another time of day and confirm the overlay still fits.

---

## 8. Reading the numbers

| Number | Where | Good value | If bad |
|---|---|---|---|
| Straightness RMS after | Lens | ≤ 0.3–0.5 px | remove non-straight lines, add border lines, change model |
| k ± σ | Lens | σ ≪ |k| | simpler model, more border lines |
| VP RMS | Vanishing pts | < 0.5 px | longer, better-fitted lines; remove outliers |
| f ± σ (MC) | Vanishing pts | σ/f < 2 % | lines with stronger convergence; use Points |
| Reprojection RMS | Points | < 0.5 px (< 1 px OK) | check coordinates, refine clicks, disable outliers |
| per-point err | Points table | green | typo in X Y Z, wrong feature |
| σ f, σ c | Points | a few px | fix parameter, add points / depth variation |
| H RMS | Plane | < 1 px | check plane points, Z = 0 only |

All σ values assume the remaining errors are random. They do **not** include systematic errors such as wrong world coordinates or the wrong distortion model — that is why visual verification matters.

---

## 9. Practise with the synthetic test scene

**Synthetic test scene** ray-traces a room corner (checkered floor 8 × 8 and two walls 8 × 5, 1-unit squares) through a camera with known parameters: 1280 × 800, `f = 880`, `c = (651.4, 402.8)`, `k1 = −0.18`, `k2 = 0.05`. Checker corners are at integer world coordinates, the origin is the corner where the floor meets both walls, X runs along the wall Y = 0, Y along the wall X = 0, Z is up.

Suggested exercise:

1. **Lens**: trace the floor–wall seams and the top wall edges (k1 k2). Expect `k1 ≈ −0.18 ± 0.01`.
2. **Vanishing points**: 3+ lines per axis along checker rows. Expect `f ≈ 880 ± 2`.
3. **Points**: 20–30 checker corners on the floor and both walls, e.g. `(3, 2, 0)`, `(0, 4, 2)`, `(5, 0, 3)`. Estimate f, pp and `k1k2`. Expect RMS ≈ 0.1 px and `f` within ~0.5 px.
4. **Verify**: the ground-truth table shows the remaining errors.

Results obtained when testing this build: Points gave `f = 879.7` (true 880), `k1 = −0.1806`, `k2 = 0.0508`, rotation error 0.04°, RMS 0.10 px. Vanishing points gave `f = 880.9 ± 1.5`. Lens gave `k1 = −0.183 ± 0.010`.

---

## 10. Keyboard & mouse reference

| Input | Action |
|---|---|
| Mouse wheel | zoom at the cursor (Ctrl + wheel / pinch = faster) |
| Right-drag, middle-drag, Space + drag | pan |
| Left-click | place a point / add a polyline vertex / start or finish a segment (tool-dependent) |
| Left-drag on empty image (Vanishing points) | draw a segment |
| Drag a handle | move it |
| ← ↑ → ↓ | nudge the selected handle 0.1 px (Shift: 1 px) |
| Delete / Backspace | delete the selected point / segment / polyline vertex |
| Enter or double-click | finish the current plumb line |
| Esc | cancel the pending segment / polyline, clear selection |
| F | fit image to window |
| ⌘Z / Ctrl-Z | undo (annotations and camera) |

---

## 11. Export format

```json
{
  "format": "manual-camera-calibration/1",
  "created": "2026-10-04T17:00:00.000Z",
  "source": "Lake Cam · 17:00:00 · via image proxy",
  "image_size": [1920, 1080],
  "pixel_convention": "OpenCV: pixel centres at integer coordinates, (0,0) = centre of top-left pixel",
  "camera_matrix": [[fx, 0, cx], [0, fy, cy], [0, 0, 1]],
  "dist_coeffs": [k1, k2, p1, p2, k3],
  "distortion_model": "opencv radial-tangential (k1,k2,p1,p2,k3)",
  "pose_valid": true,
  "rvec": [rx, ry, rz],
  "R": [[...], [...], [...]],
  "tvec": [tx, ty, tz],
  "camera_center_world": [Cx, Cy, Cz],
  "reprojection_rms_px": 0.21,
  "rms_definition": "sqrt(mean(du^2 + dv^2)) over enabled points (OpenCV convention)",
  "n_points": 18,
  "last_estimate": "points · LM",
  "std": { "f": 0.6, "cx": 0, "rot_deg": 0.03, "t": 0.01 },
  "homography_plane_to_undistorted_px": [[...], [...], [...]],
  "observations": { "points": [...], "vp_lines": [...], "plumb_lines": [...] }
}
```

`std` holds the 1-σ values of the last Points solve (`null` if none). `homography…` is `null` if Plane was not solved.

---

## 12. Using the result in OpenCV

```python
import json, numpy as np, cv2

c = json.load(open("calibration.json"))
K = np.array(c["camera_matrix"]); D = np.array(c["dist_coeffs"])
R = np.array(c["R"]); t = np.array(c["tvec"]).reshape(3, 1)
rvec = np.array(c["rvec"]).reshape(3, 1)

# project world points
X = np.array([[0, 0, 0], [1, 0, 0]], dtype=float)
uv, _ = cv2.projectPoints(X, rvec, t, K, D)

# undistort an image
img = cv2.imread("frame.jpg")
und = cv2.undistort(img, K, D)

# raw pixel -> ground plane (X, Y) with the homography
H = np.array(c["homography_plane_to_undistorted_px"])
p_und = cv2.undistortPoints(np.array([[[u, v]]], float), K, D, P=K)   # undistorted pixel
XY = cv2.perspectiveTransform(p_und, np.linalg.inv(H))
```

For Blender / three.js: vertical FOV = `2·atan(h / (2·fy))`, the principal-point offset is `((cx − (w−1)/2)/w, ((h−1)/2 − cy)/h)` in sensor fractions, and the camera-to-world rotation is `Rᵀ` combined with a flip of the y and z axes (OpenCV → OpenGL).

---

## 13. Methods (what the math does)

All in `calib.js`, plain JavaScript with no dependencies.

- **Camera model**: OpenCV pinhole + radial-tangential distortion. Undistortion is done by Newton iteration on the exact forward model (analytic 2 × 2 Jacobian, converged to 1e-14).
- **Optimisation**: Levenberg–Marquardt with Marquardt diagonal scaling and central-difference Jacobians. Covariance `σ² (JᵀJ)⁻¹` with `σ² = SSE / dof` gives the ± values.
- **Plumb lines**: minimise the perpendicular distances of undistorted points to their best-fit (total least squares) line, over `k1..k3` (tried from zero and from the current values).
- **Edge snapping**: gradient of a Gaussian-smoothed image, sampled along the edge normal, maximum with parabolic sub-pixel interpolation. Segments are fitted with TLS in undistorted coordinates.
- **Corner snapping**: iterative gradient-orthogonality solve `p = (Σ g gᵀ)⁻¹ Σ g gᵀ q`, Gaussian-weighted, rejected for edges or flat regions (eigenvalue ratio test).
- **Vanishing points**: length-weighted algebraic estimate (smallest eigenvector), then maximum-likelihood refinement on the sphere, minimising endpoint distances to the line through each segment's midpoint and the VP (handles VPs at infinity).
- **f, c from VPs**: orthogonality constraints `vᵢᵀ ω vⱼ = 0` on the image of the absolute conic. Known principal point → least-squares `f²`; three VPs → linear solve for ω (principal point = orthocentre). Rotation from `K⁻¹v`, projected onto SO(3) by polar decomposition. Uncertainty from 120 Monte-Carlo trials.
- **DLT**: Hartley-normalised 2n × 12 system → `P`, decomposed by Gram–Schmidt RQ into `K[R|t]`, with depth-sign fixing and a reflection check. With `K` known: calibrated DLT on normalised coordinates.
- **Planar case**: Hartley-normalised homography with geometric refinement. `f` from Zhang's two constraints on a single view (known principal point, square pixels), pose from `K⁻¹H`. Coplanarity is detected automatically by PCA, for any plane, not just Z = 0.
- **Final refinement**: joint LM over the parameters you selected plus `rvec` and `t`. When distortion is *kept*, it is locked in pixel space while `f` varies.

---

## 14. Troubleshooting

**"Webcam unavailable"** (your own camera)
- *Not a secure context* → serve from `http://localhost` or `https://` (see [§1](#1-files--deployment)).
- *Embedded preview not allowed* → you are inside an iframe without camera permission; open the page directly.
- *Permission denied* → you clicked Block once; re-enable it in the site settings (icon left of the address bar).
- *No camera found / in use* → connect a camera, close other apps using it.

**"Could not load …" (URL)**
- The camera is offline or the URL isn't a still image — open it in a new tab to check.
- The proxy is off and the server has no CORS or is `http://` → turn the proxy on or use your own proxy.
- The proxy can't reach the camera (private network address, IP block) → use your own proxy inside that network, or download the frame and open it as a file.

**Solve errors**
- *Inconsistent VPs / no real focal length* → a line was assigned to the wrong axis, or the axis groups aren't orthogonal in 3-D.
- *Vanishing points at infinity* → the view has no perspective along that axis; use other axes or Points.
- *Reflection / left-handed world* → swap the sign of one world axis in your coordinates.
- *Collinear points* / *underdetermined* → add points off the line, or estimate fewer parameters.
- *Points behind the camera* → wrong world coordinates or wrong axis direction.
- *Radial model folds* → too many coefficients for the data; use `k1`.

**The overlay is rotated 90° / mirrored after VPs** → **Rotate XY 90°**, or define your world axes to match.

**Annotations disappeared** → you loaded a different camera or an image of a different size (intended). Use ⌘Z or re-import your JSON.

---

## 15. Privacy

- Everything runs in your browser. No image, URL list, annotation or result is sent to any server operated by this tool.
- The only third-party request is the **image proxy fallback**, which receives the camera URL (and returns the image). Turn it off, or set the template to a proxy you control, for example a minimal Cloudflare Worker:

```js
export default { async fetch(req) {
  const u = new URL(req.url).searchParams.get('url');
  const r = await fetch(u, { cf: { cacheTtl: 0 } });
  return new Response(r.body, { headers: { 'content-type': r.headers.get('content-type') || 'image/jpeg', 'access-control-allow-origin': '*' } });
} };
```

Template: `https://your-worker.example.workers.dev/?url={url}`.

- Settings (layout, proxy, recent URLs) are kept in the browser's `localStorage`. An imported CSV list is not stored.
- If the URL list itself is private, don't deploy `webcams.csv`; import it at runtime instead.

---

## 16. Limitations

- **Single view.** Calibration uses one frame at a time. Multi-image calibration (Zhang with several checkerboard poses) would constrain intrinsics better and is the natural next step.
- **Zero skew** is assumed (true for all modern sensors). DLT reports the skew it found and then drops it.
- Plumb lines estimate radial distortion only; tangential terms (`p1, p2`) and `k3` come from Points with many well-spread points.
- Fisheye lenses beyond ~120° FOV are outside the polynomial model's comfortable range; the tool warns when the radial model folds inside the image.
- Live webcam frames must be captured before sub-pixel snapping is available.
- Accuracy is limited by your annotations and your world coordinates. The tool reports what the data supports; it cannot detect a consistently wrong measurement.
