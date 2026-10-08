import { loadPyodide } from './runtime/pyodide.mjs';

let runtime;
const status = (text) => self.postMessage({type: 'status', text});

async function initialize() {
  status('Loading the local Python runtime…');
  const pyodide = await loadPyodide({indexURL: new URL('./runtime/', import.meta.url).href});
  status('Loading numerical libraries…');
  await pyodide.loadPackage(['numpy', 'scipy', 'pyyaml']);
  const response = await fetch(new URL('./pocketseekr.zip', import.meta.url));
  if (!response.ok) throw new Error(`Unable to load PocketSeekR (${response.status})`);
  pyodide.unpackArchive(new Uint8Array(await response.arrayBuffer()), 'zip', {extractDir: '/app'});
  pyodide.runPython("import sys\nsys.path.insert(0, '/app')\nfrom pocketseekr.browser import predict_text\nimport json");
  return pyodide;
}

self.onmessage = async ({data}) => {
  if (data.type !== 'predict') return;
  try {
    runtime ??= initialize();
    const pyodide = await runtime;
    status('Searching regions and evaluating chemical support…');
    pyodide.globals.set('input_text', data.text);
    pyodide.globals.set('input_name', data.filename);
    const result = await pyodide.runPythonAsync('json.dumps(predict_text(input_text, input_name), allow_nan=False)');
    pyodide.globals.delete('input_text');
    pyodide.globals.delete('input_name');
    self.postMessage({type: 'result', id: data.id, result: JSON.parse(result)});
  } catch (error) {
    self.postMessage({type: 'error', id: data.id, message: String(error.message || error)});
  }
};
