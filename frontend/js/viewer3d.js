/**
 * 3D-визуализация расчётного блока:
 *  - геометрия блока (прозрачный параллелепипед, шкала);
 *  - скважины (ствол + интервал фильтра, цвет по назначению);
 *  - срез по выбранному слою с раскраской по полю;
 *  - линии тока с цветовой индикацией скорости.
 */
import * as THREE from 'three';
import { OrbitControls } from '../vendor/OrbitControls.js';

const COLORMAP = [
  [0.00, 0.106, 0.078, 0.392],
  [0.18, 0.122, 0.435, 0.698],
  [0.40, 0.137, 0.627, 0.541],
  [0.58, 0.353, 0.780, 0.416],
  [0.76, 0.851, 0.776, 0.180],
  [0.90, 0.910, 0.463, 0.173],
  [1.00, 0.784, 0.118, 0.118],
];

export function colorFor(t) {
  t = Math.max(0, Math.min(1, t));
  const n = COLORMAP.length - 1;
  const i = Math.min(Math.floor(t * n), n - 1);
  const f = t * n - i;
  const a = COLORMAP[i], b = COLORMAP[i + 1];
  return [
    (a[0] + (b[0] - a[0]) * f) * 255,
    (a[1] + (b[1] - a[1]) * f) * 255,
    (a[2] + (b[2] - a[2]) * f) * 255,
  ];
}

export function cssColor(t) {
  const c = colorFor(t);
  return `rgb(${c[0] | 0},${c[1] | 0},${c[2] | 0})`;
}

export function webglAvailable() {
  try {
    const canvas = document.createElement('canvas');
    return !!(window.WebGLRenderingContext
      && (canvas.getContext('webgl2') || canvas.getContext('webgl')));
  } catch (e) {
    return false;
  }
}

export class Viewer3D {
  constructor(container) {
    this.container = container;
    this.verticalExaggeration = 4.0;
    this.available = false;

    if (!webglAvailable()) {
      container.innerHTML = `
        <div style="display:flex;align-items:center;justify-content:center;height:100%;
                    text-align:center;padding:24px;color:#7d8b9a;font-size:13px;line-height:1.6">
          <div>
            <div style="font-size:30px;margin-bottom:10px">◈</div>
            <b style="color:#d6dee8">3D-вид недоступен</b><br>
            В браузере отключён WebGL.<br>
            Включите аппаратное ускорение или используйте карту среза
            и кривые в правой панели — расчёты работают полностью.
          </div>
        </div>`;
      return;
    }

    this.blockGroup = new THREE.Group();
    this.wellGroup = new THREE.Group();
    this.sliceGroup = new THREE.Group();
    this.lineGroup = new THREE.Group();

    this.scene = new THREE.Scene();
    this.scene.background = new THREE.Color(0x0e1319);
    this.scene.add(this.blockGroup, this.wellGroup, this.sliceGroup, this.lineGroup);

    const w = container.clientWidth || 800;
    const h = container.clientHeight || 600;
    this.camera = new THREE.PerspectiveCamera(45, w / h, 1, 20000);
    this.camera.position.set(700, -700, 500);
    this.camera.up.set(0, 0, 1);

    this.renderer = new THREE.WebGLRenderer({ antialias: true });
    this.renderer.setSize(w, h);
    this.renderer.setPixelRatio(Math.min(devicePixelRatio, 2));
    container.appendChild(this.renderer.domElement);

    this.available = true;
    this.controls = new OrbitControls(this.camera, this.renderer.domElement);
    this.controls.enableDamping = true;
    this.controls.dampingFactor = 0.08;

    this.scene.add(new THREE.AmbientLight(0xffffff, 0.75));
    const dir = new THREE.DirectionalLight(0xffffff, 0.9);
    dir.position.set(300, -400, 600);
    this.scene.add(dir);
    const dir2 = new THREE.DirectionalLight(0xffffff, 0.4);
    dir2.position.set(-400, 300, 200);
    this.scene.add(dir2);

    this.gridHelper = null;
    this._animate();
  }

  resize() {
    if (!this.available) return;
    const w = this.container.clientWidth;
    const h = this.container.clientHeight;
    if (!w || !h) return;
    this.camera.aspect = w / h;
    this.camera.updateProjectionMatrix();
    this.renderer.setSize(w, h);
  }

  _animate() {
    if (!this.available) return;
    requestAnimationFrame(() => this._animate());
    this.controls.update();
    this.renderer.render(this.scene, this.camera);
  }

  clear(group) {
    if (!group) return;
    while (group.children.length) {
      const c = group.children.pop();
      c.traverse?.((o) => {
        if (o.geometry) o.geometry.dispose();
        if (o.material) {
          if (Array.isArray(o.material)) o.material.forEach((m) => m.dispose());
          else o.material.dispose();
        }
      });
    }
  }

  setBlock(grid) {
    if (!this.available) { this.grid = grid; return; }
    this.clear(this.blockGroup);
    this.grid = grid;
    const { nx, ny, nz, dx, dy, dz } = grid;
    const w = nx * dx, h = ny * dy, d = nz * dz;
    const geo = new THREE.BoxGeometry(w, h, d * this.verticalExaggeration);
    const edges = new THREE.EdgesGeometry(geo);
    const line = new THREE.LineSegments(edges, new THREE.LineBasicMaterial({
      color: 0x5c6b7a, transparent: true, opacity: 0.9,
    }));
    line.position.set(
      grid.x0 + w / 2,
      grid.y0 + h / 2,
      grid.z0 + (d * this.verticalExaggeration) / 2,
    );
    this.blockGroup.add(line);

    const face = new THREE.Mesh(
      geo,
      new THREE.MeshBasicMaterial({
        color: 0x2b3a4a, transparent: true, opacity: 0.06, side: THREE.DoubleSide,
        depthWrite: false,
      }),
    );
    face.position.copy(line.position);
    this.blockGroup.add(face);

    // Ось X на поверхности (м)
    this._axisLine(new THREE.Vector3(grid.x0, grid.y0, grid.z0),
      new THREE.Vector3(grid.x0 + w, grid.y0, grid.z0), 0xd44);
    this._axisLine(new THREE.Vector3(grid.x0, grid.y0, grid.z0),
      new THREE.Vector3(grid.x0, grid.y0 + h, grid.z0), 0x4d8);
    this._axisLine(new THREE.Vector3(grid.x0, grid.y0, grid.z0),
      new THREE.Vector3(grid.x0, grid.y0, grid.z0 + d * this.verticalExaggeration), 0x8d4);

    this.fitCamera();
  }

  _axisLine(a, b, color) {
    const g = new THREE.BufferGeometry().setFromPoints([a, b]);
    this.blockGroup.add(new THREE.Line(g, new THREE.LineBasicMaterial({ color })));
  }

  fitCamera() {
    if (!this.available || !this.grid) return;
    const { nx, ny, nz, dx, dy, dz, x0, y0, z0 } = this.grid;
    const w = nx * dx, h = ny * dy, d = nz * dz * this.verticalExaggeration;
    const cx = x0 + w / 2;
    const cy = y0 + h / 2;
    const cz = z0 + d / 2;
    const radius = Math.sqrt(w * w + h * h + d * d) / 2;
    const fov = (this.camera.fov * Math.PI) / 180;
    const dist = (radius / Math.sin(fov / 2)) * 1.05;
    const dir = new THREE.Vector3(0.62, -0.72, 0.52).normalize();
    this.camera.position.copy(this.controls.target.clone().set(cx, cy, cz)).add(dir.multiplyScalar(dist));
    this.controls.target.set(cx, cy, cz);
    this.controls.update();
  }

  setWells(wells, grid) {
    if (!this.available) return;
    this.clear(this.wellGroup);
    if (!grid) return;
    const ve = this.verticalExaggeration;
    const radius = Math.max(1.2, Math.min(grid.dx, grid.dy) * 0.045);
    const colorOf = (t) => (t === 'injection' ? 0x3d8bff : t === 'production' ? 0xff8c3a : 0x9aa7b4);

    for (const w of wells) {
      const col = colorOf(w.well_type);
      const x = w.x, y = w.y;
      const zHead = (w.z ?? grid.z0) * 1 + grid.z0;
      const top = new THREE.Vector3(x, y, grid.z0);
      const bottom = new THREE.Vector3(x, y, grid.z0 + 0.5 * grid.dz * ve);

      const stemH = Math.max(Math.abs(top.z - bottom.z), 0.5);
      const stem = new THREE.Mesh(
        new THREE.CylinderGeometry(radius * 0.6, radius * 0.6, stemH, 8),
        new THREE.MeshBasicMaterial({ color: col }),
      );
      stem.rotation.x = Math.PI / 2;
      stem.position.set(x, y, Math.min(top.z, bottom.z) + stemH / 2);
      this.wellGroup.add(stem);

      const screens = [];
      if (w.screen_top != null && w.screen_bottom != null) {
        screens.push([w.screen_top, w.screen_bottom]);
      }
      (w.screens || []).forEach((s) => screens.push([s.top, s.bottom]));
      for (const [topD, botD] of screens) {
        const z1 = grid.z0 + (topD - grid.z0) * ve;
        const z2 = grid.z0 + (botD - grid.z0) * ve;
        const h = Math.max(z2 - z1, 1);
        const scr = new THREE.Mesh(
          new THREE.CylinderGeometry(radius * 1.5, radius * 1.5, h, 10),
          new THREE.MeshBasicMaterial({ color: col, transparent: true, opacity: 0.55 }),
        );
        scr.rotation.x = Math.PI / 2;
        scr.position.set(x, y, z1 + h / 2);
        this.wellGroup.add(scr);
      }

      // Указатель на устье
      const marker = new THREE.Mesh(
        new THREE.ConeGeometry(radius * 2.2, radius * 3.2, 10),
        new THREE.MeshBasicMaterial({ color: col }),
      );
      marker.position.set(x, y, grid.z0 + radius * 2.0);
      this.wellGroup.add(marker);
    }
  }

  setSlice(data, nx, ny, grid, layer, opts = {}) {
    if (!this.available) return;
    this.clear(this.sliceGroup);
    if (!data || !grid) return;
    const ve = this.verticalExaggeration;
    const w = grid.nx * grid.dx;
    const h = grid.ny * grid.dy;

    // Раскраска поля в RGBA (цветовая шкала)
    const min = opts.min ?? 0;
    const max = opts.max ?? 1;
    const span = (max - min) || 1;
    const rgba = new Uint8Array(nx * ny * 4);
    for (let j = 0; j < ny; j++) {
      for (let i = 0; i < nx; i++) {
        const v = data[i * ny + j];
        const t = Math.max(0, Math.min(1, (v - min) / span));
        const c = colorFor(t);
        // строка 0 текстуры соответствует югу (меньшее Y)
        const row = j;
        const idx = (row * nx + i) * 4;
        rgba[idx] = c[0]; rgba[idx + 1] = c[1]; rgba[idx + 2] = c[2]; rgba[idx + 3] = 255;
      }
    }

    const tex = new THREE.DataTexture(rgba, nx, ny, THREE.RGBAFormat);
    tex.needsUpdate = true;
    tex.minFilter = THREE.LinearFilter;
    tex.magFilter = THREE.LinearFilter;
    tex.colorSpace = THREE.SRGBColorSpace;
    tex.flipY = false;

    const mat = new THREE.MeshBasicMaterial({
      map: tex, transparent: true, opacity: opts.opacity ?? 0.85,
      side: THREE.DoubleSide, depthWrite: false,
    });
    // Срез — горизонтальная плоскость постоянной глубины (в системе координат
    // ось Z направлена вниз, поэтому плоскость XY не поворачиваем)
    const plane = new THREE.Mesh(new THREE.PlaneGeometry(w, h), mat);
    const z = grid.z0 + (layer + 0.5) * grid.dz * ve;
    plane.position.set(grid.x0 + w / 2, grid.y0 + h / 2, z);
    this.sliceGroup.add(plane);
    this.sliceTexture = tex;
  }

  setStreamlines(lines, grid) {
    if (!this.available) return;
    this.clear(this.lineGroup);
    if (!lines || !grid) return;
    const ve = this.verticalExaggeration;
    const positions = [];
    const colors = [];
    let vmax = 1e-9;
    for (const l of lines) {
      for (const s of l.speeds || []) vmax = Math.max(vmax, s);
    }
    const stride = 2;
    for (const l of lines) {
      const pts = l.points || [];
      for (let i = 0; i + stride < pts.length; i += stride) {
        const a = pts[i], b = pts[i + stride];
        positions.push(a[0], a[1], grid.z0 + (a[2] - grid.z0) * ve);
        positions.push(b[0], b[1], grid.z0 + (b[2] - grid.z0) * ve);
        const sp = l.speeds || [];
        const t = Math.min(1, (sp[i] || 0) / vmax);
        const c = colorFor(t);
        colors.push(c[0] / 255, c[1] / 255, c[2] / 255);
        colors.push(c[0] / 255, c[1] / 255, c[2] / 255);
      }
    }
    if (!positions.length) return;
    const geo = new THREE.BufferGeometry();
    geo.setAttribute('position', new THREE.Float32BufferAttribute(positions, 3));
    geo.setAttribute('color', new THREE.Float32BufferAttribute(colors, 3));
    const mat = new THREE.LineBasicMaterial({
      vertexColors: true, transparent: true, opacity: 0.85,
    });
    this.lineGroup.add(new THREE.LineSegments(geo, mat));
    this.lineVmax = vmax;
  }

  setVisible(what, value) {
    if (!this.available) return;
    if (what === 'wells') this.wellGroup.visible = value;
    if (what === 'lines') this.lineGroup.visible = value;
    if (what === 'slice') this.sliceGroup.visible = value;
  }

  setSliceOpacity(v) {
    if (!this.available) return;
    this.sliceGroup.traverse((o) => {
      if (o.material && o.material.map) o.material.opacity = v;
    });
  }
}
