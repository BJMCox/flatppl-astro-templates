// Resolve dependencies from disk and evaluate through the JS module loader.
const fs = require('node:fs');
const path = require('node:path');
const root = path.resolve(__dirname, '..');
const engine = require(path.resolve(process.env.FLATPPL_JS_DIR || path.join(root, '../flatppl-js'), 'packages/engine/index.ts'));
const request = JSON.parse(fs.readFileSync(0, 'utf8'));
const expressions = Array.isArray(request) ? request : null;
const source = expressions ? 'flatppl_compat = "0.1"\na = load_module("../main.flatppl")\n'
  + expressions.map((expression, i) => `probe${i} = ${expression}`).join('\n') : request.source;
const primary = path.join(__dirname, 'query.flatppl');
const sources = {};
function walk(importer, text) {
  for (const dep of engine.moduleDeps(text)) {
    const resolved = engine.moduleResolve.resolveModulePath(importer, dep);
    if (sources[resolved]) continue;
    let local = resolved;
    const prefix = 'https://raw.githubusercontent.com/flatppl-community/flatppl-stats-templates/';
    if (resolved.startsWith(prefix) && process.env.FLATPPL_STATS_TEMPLATES_DIR) {
      const relative = resolved.slice(prefix.length).split('/').slice(1).join('/');
      local = path.join(process.env.FLATPPL_STATS_TEMPLATES_DIR, relative);
    }
    sources[resolved] = fs.readFileSync(local, 'utf8');
    walk(resolved, sources[resolved]);
  }
}
walk(primary, source);
const model = engine.processSource(source, { path: primary, bundle: { sources } });
const errors = model.diagnostics.filter(d => d.severity === 'error');
if (errors.length) throw new Error(JSON.stringify(errors));
const result = engine.orchestrator.buildDerivations(model.linkedBindings, { moduleRegistry: model.linkedModuleRegistry });
const { createWorkerHandler } = require(path.resolve(process.env.FLATPPL_JS_DIR || path.join(root, '../flatppl-js'), 'packages/engine/worker.ts'));
const worker = createWorkerHandler();
worker.handle({ type: 'init', seed: 1 });
const cache = new Map();
const context = {
  ...result,
  sampleCount: 1,
  rootSeed: 1,
  getMeasure(name) {
    if (!cache.has(name)) cache.set(name, engine.materialiser.materialiseMeasure(name, context));
    return cache.get(name);
  },
  async sendWorker(message) {
    const reply = await worker.handle(message);
    if (reply?.type === 'error') throw new Error(reply.message);
    return reply;
  },
};
const names = expressions ? expressions.map((_, i) => `probe${i}`) : request.bindings;
function plain(value) {
  if (Array.isArray(value)) return value.map(plain);
  if (value && value.shape && value.data) {
    return value.shape.length ? Array.from(value.data, plain) : plain(value.data[0]);
  }
  return value;
}
Promise.all(names.map(async name => {
  if (result.fixedValues.has(name)) return plain(result.fixedValues.get(name));
  const measure = await context.getMeasure(name);
  if (!measure?.samples) throw new Error(`${name} did not evaluate`);
  return plain(measure.samples[0]);
})).then(values => process.stdout.write(JSON.stringify(values))).catch(error => {
  process.stderr.write(String(error));
  process.exitCode = 1;
});
