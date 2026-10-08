const $ = (id) => document.getElementById(id);
const colors = ['#21867a', '#cd765d', '#9574ad', '#bd953e', '#548ab3'];
let inputText = '', inputName = '', worker, report, selected = 0, viewer;
let busy = false, loadingInput = false, requestId = 0, started = 0, timer;

function setStatus(message, kind = '') {
  $('status').textContent = message;
  $('status-box').className = `status-box ${kind}`;
}

function clearResults() {
  report = undefined;
  $('results').hidden = true;
  $('placeholder').hidden = false;
  $('reset-view').disabled = true;
  $('viewer-title').textContent = 'RNA structure viewer';
  $('viewer-info').textContent = 'Drag to rotate · scroll to zoom';
  if (viewer) { viewer.clear(); viewer.render(); }
}

function setInput(text, name) {
  inputText = text;
  inputName = name;
  clearResults();
  $('filename').textContent = name;
  $('viewer-title').textContent = 'RNA structure viewer';
  $('run').disabled = !text || loadingInput;
  $('elapsed').textContent = '';
  setStatus('Structure selected. Ready to identify pockets.');
}

function lockInput(loading) {
  loadingInput = loading;
  $('structure').disabled = busy || loading;
  document.querySelectorAll('[data-example]').forEach(button => { button.disabled = busy || loading; });
  $('run').disabled = busy || loading || !inputText;
}

async function readFile(file) {
  if (!file || busy || loadingInput) return;
  clearResults();
  inputText = '';
  $('filename').textContent = 'No structure selected';
  $('run').disabled = true;
  if (!file.name.toLowerCase().endsWith('.pdb')) {
    setStatus('Please choose a PDB file. Direct mmCIF input is not supported.', 'error');
    return;
  }
  if (file.size > 20 * 1024 * 1024) {
    setStatus('The browser accepts files up to 20 MB. Use the Python CLI for larger files.', 'error');
    return;
  }
  lockInput(true);
  try { setInput(await file.text(), file.name); }
  catch (error) { setStatus(`Unable to read the file: ${error.message}`, 'error'); }
  finally { lockInput(false); }
}

function finish() {
  busy = false;
  clearInterval(timer);
  $('run').disabled = !inputText;
  $('structure').disabled = false;
  document.querySelectorAll('[data-example]').forEach(button => { button.disabled = false; });
  $('cancel').hidden = true;
}

function download(content, name, type) {
  const url = URL.createObjectURL(new Blob([content], {type}));
  const link = document.createElement('a');
  link.href = url;
  link.download = name;
  link.click();
  setTimeout(() => URL.revokeObjectURL(url), 1000);
}

function stem() { return inputName.replace(/\.pdb$/i, '').replace(/[^a-zA-Z0-9_.-]/g, '_') || 'rna'; }

function showPocket(index) {
  selected = index;
  const pocket = report.pockets[index];
  viewer.setStyle({}, {cartoon: {color: '#aab6c3', opacity: .7}, stick: {radius: .08, color: '#b8c4d0', opacity: .4}});
  viewer.setStyle({serial: pocket.atom_indices.map((i) => i + 1)}, {stick: {radius: .18, color: colors[index]}, cartoon: {color: colors[index]}});
  viewer.render();
  document.querySelectorAll('.pocket-card').forEach((button, i) => {
    button.classList.toggle('active', i === index);
    button.setAttribute('aria-pressed', String(i === index));
    button.closest('tr').classList.toggle('selected', i === index);
  });
  const residues = new Set(pocket.crop_atom_identities.map((a) => JSON.stringify([a.chain_id, a.residue_id, a.residue_name]))).size;
  $('crop-info').textContent = `Pocket ${pocket.rank} · ${residues} observed nucleotides · ${pocket.atom_indices.length} RNA atoms`;
  $('download-crop').disabled = false;
}

function showResults(result) {
  report = result;
  $('placeholder').hidden = true;
  $('viewer-title').textContent = inputName;
  $('viewer-info').textContent = `${result.rna_atom_count.toLocaleString()} RNA atoms · rotate by dragging`;
  if (!globalThis.$3Dmol) throw new Error('The 3D viewer could not be loaded.');
  viewer ??= globalThis.$3Dmol.createViewer($('viewer'), {backgroundColor: '#ffffff'});
  viewer.clear();
  viewer.addModel(result.rna_pdb, 'pdb');
  viewer.setStyle({}, {cartoon: {color: '#aab6c3'}, stick: {radius: .1, color: '#aab6c3'}});
  for (let i = 0; i < result.pockets.length; i++) {
    const [x, y, z] = result.pockets[i].center;
    viewer.addSphere({center: {x, y, z}, radius: .85, color: colors[i]});
    viewer.addLabel(String(i + 1), {position: {x: x + 1.2, y: y + 1.2, z}, fontSize: 14, fontColor: colors[i], backgroundOpacity: 0});
  }
  viewer.zoomTo();
  viewer.render();
  $('reset-view').disabled = false;
  $('pocket-list').replaceChildren();
  $('pocket-count').textContent = `${result.pockets.length} / 5`;
  $('results').hidden = false;
  $('download-crop').disabled = true;
  $('crop-info').textContent = '';
  for (let i = 0; i < result.pockets.length; i++) {
    const pocket = result.pockets[i];
    const card = document.createElement('button');
    card.type = 'button';
    card.className = 'pocket-card';
    const dot = document.createElement('span');
    dot.className = 'pocket-dot'; dot.style.background = colors[i];
    card.append(dot, document.createTextNode(`Pocket ${pocket.rank}`));
    card.addEventListener('click', () => showPocket(i));
    const row = document.createElement('tr');
    const choice = document.createElement('td');
    choice.append(card); row.append(choice);
    for (const value of [pocket.score.toFixed(3), pocket.selection_gain.toFixed(3),
      pocket.center.map((v) => v.toFixed(2)).join(', '), String(pocket.atom_indices.length)]) {
      const cell = document.createElement('td'); cell.textContent = value; row.append(cell);
    }
    $('pocket-list').append(row);
  }
  if (result.pockets.length) showPocket(0);
  else $('crop-info').textContent = 'No eligible pocket was found. The method does not create replacement candidates.';
  window.dispatchEvent(new Event('resize'));
}

function startWorker() {
  worker = new Worker(new URL('./worker.js', import.meta.url), {type: 'module'});
  const activeWorker = worker;
  worker.onmessage = ({data}) => {
    if (worker !== activeWorker) return;
    if (data.type === 'status' && busy) setStatus(data.text, 'running');
    if (data.id !== requestId) return;
    if (data.type === 'result') {
      try { showResults(data.result); setStatus(`Completed · ${data.result.pockets.length} pockets identified.`); }
      catch (error) { clearResults(); setStatus(error.message, 'error'); }
      finish();
    } else if (data.type === 'error') {
      clearResults(); setStatus(data.message, 'error');
      worker.terminate(); worker = undefined;
      finish();
    }
  };
  worker.onerror = (event) => {
    if (worker !== activeWorker) return;
    clearResults(); setStatus(`Runtime error: ${event.message || 'Unable to start browser computation.'}`, 'error');
    worker.terminate(); worker = undefined;
    finish();
  };
}

$('structure').addEventListener('change', (event) => readFile(event.target.files[0]));
const zone = document.querySelector('.file-zone');
zone.addEventListener('dragover', (event) => { event.preventDefault(); if (!busy && !loadingInput) zone.classList.add('dragging'); });
zone.addEventListener('dragleave', () => zone.classList.remove('dragging'));
zone.addEventListener('drop', (event) => { event.preventDefault(); zone.classList.remove('dragging'); readFile(event.dataTransfer.files[0]); });
async function loadExample() {
  if (busy || loadingInput) return;
  lockInput(true);
  try {
    const response = await fetch('./examples/1F1T_rna.pdb');
    if (!response.ok) throw new Error(`Example download failed (${response.status})`);
    setInput(await response.text(), '1F1T_rna.pdb');
    $('submit').scrollIntoView({block: 'start'});
  } catch (error) { setStatus(error.message, 'error'); }
  finally { lockInput(false); }
}
document.querySelectorAll('[data-example]').forEach(button => button.addEventListener('click', loadExample));
$('run').addEventListener('click', () => {
  if (busy || loadingInput || !inputText) return;
  busy = true; clearResults(); requestId++;
  lockInput(false); $('cancel').hidden = false;
  started = performance.now();
  const update = () => { $('elapsed').textContent = `Total elapsed: ${Math.floor((performance.now() - started) / 1000)} s`; };
  update(); timer = setInterval(update, 1000);
  setStatus('Starting local computation…', 'running');
  if (!worker) startWorker();
  worker.postMessage({type: 'predict', id: requestId, text: inputText, filename: inputName});
});
$('cancel').addEventListener('click', () => {
  requestId++; worker.terminate(); worker = undefined;
  finish(); setStatus('Calculation cancelled. You can run again.');
});
$('reset-view').addEventListener('click', () => { viewer.zoomTo(); viewer.render(); });
$('download-json').addEventListener('click', () => download(JSON.stringify(report, null, 2) + '\n', `${stem()}_pockets.json`, 'application/json'));
$('download-rna').addEventListener('click', () => download(report.rna_pdb, `${stem()}_rna.pdb`, 'chemical/x-pdb'));
$('download-crop').addEventListener('click', () => download(report.pockets[selected].crop_pdb, `${stem()}_pocket_${selected + 1}.pdb`, 'chemical/x-pdb'));
window.addEventListener('resize', () => { if (viewer) viewer.resize(); });
