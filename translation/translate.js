const OpenApi = require('@alicloud/openapi-client');
const Util = require('@alicloud/tea-util');
const commander = require('commander');
const fs = require('fs');
const { exit } = require('process');

const dictionaryFile = 'dictionary.txt';

const program = new commander.Command();

program
  .version('1.0.0', '-v, --version')
  .usage('[OPTIONS]...')
  .requiredOption('-f, --input <file-path>', 'Input file path')
  .requiredOption('-s, --secret <value>', 'OSS Access key secret')
  .requiredOption('-i, --id <value>', 'OSS Access key ID')
  .requiredOption('-e, --endpoint <value>', 'MT endpoint')
  .requiredOption('-o, --output <file-path>', 'Output file path')
  .parse();

const options = program.opts();

const accessKeyId = options.id;
const accessKeySecret = options.secret;
const endpoint = options.endpoint;

let translatedDictionary = [];

const sourceLanguage = 'en';
const targetLanguage = 'zh';

const config = new OpenApi.Config({
  accessKeyId: accessKeyId,
  accessKeySecret: accessKeySecret,
  endpoint: endpoint,
});

function loadDictionary() {
  const fileData = fs
    .readFileSync(`${dictionaryFile}`, 'utf-8')
    .replace(/\r|\r\n|\n/g, '\n')
    .split('\n');

  return fileData
    .map((line) => {
      const [source, translated] = line.split('\t');
      return { source, translated };
    })
    .filter((datum) => datum.source && datum.translated);
}

const client = new OpenApi.default(config);

async function run() {
  const data = loadData();
  translatedDictionary = loadDictionary();

  // console.log(translatedDictionary);

  const total = data.length;
  console.log('Total:', total);

  // Get translate data
  for (const index in data) {
    const datum = data[index];

    const translatedInfo = await translate(datum['source']);
    datum['result'] = translatedInfo;

    // Write result to file
    writeResultToFile(datum);
  }
}

run();

function writeResultToFile(datum) {
  const headers = ['source', 'result'];
  const tsvData = headers.map((header) => datum[header]).join('\t');
  fs.appendFileSync(options.output, `${tsvData}\n`, 'utf8');
}

function loadData() {
  const headers = ['source'];
  const fileData = fs
    .readFileSync(options.input, 'utf-8')
    .replace(/\r|\r\n|\n/g, '\n')
    .split('\n')
    .filter((d) => d);

  return fileData.map((element) => {
    const datum = {};
    element.split('\t').forEach((item, index) => {
      datum[headers[index]] = item;
    });
    return datum;
  });
}

function writeTranslatedToDictionary(source, translated) {
  // Check if the source text is already existed in the dictionary
  const translatedResult = findInDictionary(source);
  if (translatedResult) {
    return;
  }

  const sourceText = source.toLowerCase();
  translatedDictionary.push({
    source: sourceText,
    translated: translated,
  });
  fs.appendFileSync(
    `${dictionaryFile}`,
    `${sourceText}\t${translated}\n`,
    'utf8'
  );
}

function findInDictionary(source) {
  const result = translatedDictionary.find(
    (d) => d.source.toLowerCase() == source.toLowerCase()
  );
  if (result) {
    return result.translated;
  }
  return null;
}

function createApiInfo(action) {
  const params = new OpenApi.Params({
    // 接口名称
    action,
    // 接口版本
    version: '2018-10-12',
    // 接口协议
    protocol: 'HTTPS',
    // 接口 HTTP 方法
    method: 'POST',
    authType: 'AK',
    style: 'RPC',
    // 接口 PATH
    pathname: `/`,
    // 接口请求体内容格式
    reqBodyType: 'formData',
    // 接口响应体内容格式
    bodyType: 'json',
  });
  return params;
}

async function wait(ms) {
  return new Promise((resolve) => setTimeout(resolve, ms));
}

async function translate(sourceText) {
  const inputSource = sourceText.toLowerCase();

  const dictTranslated = findInDictionary(inputSource);

  if (dictTranslated) {
    return dictTranslated;
  }

  const params = createApiInfo('TranslateGeneral');
  const body = {
    FormatType: 'text',
    SourceLanguage: sourceLanguage,
    TargetLanguage: targetLanguage,
    SourceText: inputSource,
    Scene: 'general', //  'general' || 'medical'
  };

  let result;

  try {
    const request = new OpenApi.OpenApiRequest({
      body: body,
    });
    const runtime = new Util.RuntimeOptions({});
    // Wait 3 seconds before perform the next api call
    await wait(2000);
    result = await client.callApi(params, request, runtime);
    const translated = result.body.Data.Translated;
    writeTranslatedToDictionary(inputSource, translated);
    return translated;
  } catch (error) {
    console.error('Error calling Aliyun Translate API:', error.message);
    console.error(result ? result.body : 'NO RESULT');
  }
}
