/* Exercise the actual WASM core through the public UI, also at a Pages subpath. */
const assert = require('node:assert/strict');
const fs = require('node:fs/promises');
const path = require('node:path');
const http = require('node:http');
const {execFileSync} = require('node:child_process');
const {chromium} = require('playwright');

const root = path.resolve(__dirname, '..');
const output = path.join(root, 'outputs/browser-verification');
const site = path.join(root, 'outputs/site');
const mime = {'.html':'text/html', '.js':'text/javascript', '.mjs':'text/javascript', '.css':'text/css', '.wasm':'application/wasm'};
let maximumError = 0;
function compare(actual, expected, key = '') {
  if (typeof expected === 'number') {
    assert.equal(typeof actual, 'number', key);
    if (Number.isInteger(expected)) assert.equal(actual, expected, key);
    else {
      const error = Math.abs(actual - expected);
      maximumError = Math.max(maximumError, error);
      assert.ok(error <= 1e-9 * Math.max(1, Math.abs(expected)), `${key}: ${actual} != ${expected}`);
    }
  } else if (expected && typeof expected === 'object') {
    assert.deepEqual(Object.keys(actual).sort(), Object.keys(expected).sort(), key);
    for (const k of Object.keys(expected)) compare(actual[k], expected[k], `${key}.${k}`);
  } else assert.equal(actual, expected, key);
}

(async () => {
  await fs.mkdir(output, {recursive:true});
  let server, browser;
  const pageErrors = [], requests = [], verification = [];
  try {
    let url = process.env.SITE_URL;
    if (!url) {
      server = http.createServer(async (request, response) => {
        try {
          assert.equal(request.method, 'GET');
          const name = decodeURIComponent(new URL(request.url, 'http://localhost').pathname);
          assert.ok(name.startsWith('/PocketSeekR/'));
          const file = path.resolve(site, name.slice('/PocketSeekR/'.length) || 'index.html');
          assert.ok(file.startsWith(site + path.sep));
          response.setHeader('Content-Type', mime[path.extname(file)] || 'application/octet-stream');
          response.end(await fs.readFile(file));
        } catch { response.writeHead(404); response.end(); }
      });
      await new Promise(resolve => server.listen(0, '127.0.0.1', resolve));
      url = `http://127.0.0.1:${server.address().port}/PocketSeekR/`;
    }
    const origin = new URL(url).origin;
    browser = await chromium.launch({headless:true,
      ...(process.env.CHROME_PATH ? {executablePath:process.env.CHROME_PATH} : {}),
      args:['--use-gl=angle', '--use-angle=swiftshader', '--enable-unsafe-swiftshader']});
    const page = await browser.newPage({viewport:{width:1440,height:1100}, acceptDownloads:true});
    page.on('pageerror', error => pageErrors.push(error.message));
    page.on('request', r => requests.push({url:r.url(),method:r.method()}));
    await page.goto(url);
    await page.locator('#example').click();
    await page.waitForFunction(() => !document.getElementById('run').disabled);
    // A killed worker must not interfere with a subsequent calculation.
    await page.locator('#run').click();
    await page.locator('#cancel').click();
    assert.match(await page.locator('#status').textContent(), /cancelled/);
    assert.equal(await page.locator('#run').isEnabled(), true);

    const inputs = ['1F1T', '1NTA'];
    for (const entry of inputs) {
      const source = path.join(root, 'tests/data', entry+'.pdb');
      await page.locator('#structure').setInputFiles(source);
      await page.waitForFunction(() => !document.getElementById('run').disabled);
      const started = Date.now();
      await page.locator('#run').click();
      await page.waitForFunction(() => document.getElementById('cancel').hidden, {timeout:180000});
      assert.match(await page.locator('#status').textContent(), /^Completed/, `${entry}: ${await page.locator('#status').textContent()}`);
      const downloadEvent = page.waitForEvent('download');
      await page.locator('#download-json').click();
      const download = await downloadEvent;
      const destination = path.join(output, entry+'.json');
      await download.saveAs(destination);
      const actual = JSON.parse(await fs.readFile(destination, 'utf8'));
      const expected = JSON.parse(execFileSync(process.env.PYTHON || 'python3', ['-c',
        'import json,sys; from pathlib import Path; from pocketseekr.browser import predict_text; print(json.dumps(predict_text(Path(sys.argv[1]).read_text(),Path(sys.argv[1]).name)))', source],
        {cwd:root, env:{...process.env,PYTHONPATH:path.join(root,'src')}, maxBuffer:30*1024*1024}).toString());
      compare(actual, expected, entry);
      assert.equal(await page.locator('.pocket-card').count(), 5);
      assert.ok(await page.locator('#viewer canvas').count() > 0);
      await page.locator('.pocket-card').nth(2).click();
      assert.match(await page.locator('#crop-info').textContent(), /Pocket 3/);
      const cropEvent = page.waitForEvent('download');
      await page.locator('#download-crop').click();
      const crop = await cropEvent; await crop.saveAs(path.join(output, entry+'_crop.pdb'));
      assert.equal(await fs.readFile(path.join(output, entry+'_crop.pdb'), 'utf8'), actual.pockets[2].crop_pdb);
      verification.push({pdb:entry, regions:actual.region_pool.length, pockets:actual.pockets.length,
        full_report_matches_native_with_tolerance:true, seconds:(Date.now()-started)/1000});
      if (entry === '1F1T') await page.screenshot({path:path.join(output,'desktop.png'),fullPage:true});
    }
    await page.setViewportSize({width:390,height:844});
    await page.screenshot({path:path.join(output,'mobile.png'),fullPage:true});
    assert.equal(await page.evaluate(() => document.documentElement.scrollWidth > innerWidth), false);
    await page.locator('#structure').setInputFiles({name:'empty.pdb',mimeType:'chemical/x-pdb',buffer:Buffer.from('END\n')});
    await page.waitForFunction(() => !document.getElementById('run').disabled);
    await page.locator('#run').click();
    await page.waitForFunction(() => document.getElementById('cancel').hidden, {timeout:180000});
    assert.match(await page.locator('#status').textContent(), /no RNA atoms/);
    assert.equal(await page.locator('#results').isVisible(), false);
    // Errors restart the runtime cleanly; one eligible pocket is not padded to five.
    await page.locator('#structure').setInputFiles(path.join(root,'examples/synthetic_shell.pdb'));
    await page.waitForFunction(() => !document.getElementById('run').disabled);
    await page.locator('#run').click();
    await page.waitForFunction(() => document.getElementById('cancel').hidden, {timeout:180000});
    assert.match(await page.locator('#status').textContent(), /^Completed/);
    assert.equal(await page.locator('.pocket-card').count(), 1);
    assert.deepEqual(pageErrors, []);
    assert.ok(requests.every(r => new URL(r.url).origin === origin && r.method === 'GET'), 'Structure data must never be sent to a service');
    const record = {browser:await browser.version(), url, comparison_tolerance:1e-9, maximum_absolute_numeric_error:maximumError,
      comparisons:verification, checks:['Pages subpath','cancel and rerun','real PDB upload','JSON export','crop export','3D viewer','mobile overflow','invalid RNA error','recovery after error','no candidate padding','all requests same-origin GET']};
    await fs.writeFile(path.join(output,'summary.json'), JSON.stringify(record,null,2)+'\n');
    console.log(JSON.stringify(record,null,2));
  } finally {
    if (browser) await browser.close();
    if (server) await new Promise(resolve => server.close(resolve));
  }
})().catch(error => {console.error(error);process.exitCode=1;});
