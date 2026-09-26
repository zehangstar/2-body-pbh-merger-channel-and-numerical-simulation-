"""Build a versioned project archive and a Colab launcher from this checkout."""
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import subprocess
import zipfile

ROOT = Path(__file__).resolve().parents[1]
BUILD = ROOT / 'cloud' / 'build'
DRIVE_NAME = 'PBH_Merger_Colab'


def digest(data):
    return hashlib.sha256(data).hexdigest()


def main():
    BUILD.mkdir(parents=True, exist_ok=True)
    revision = subprocess.check_output(['git', 'rev-parse', 'HEAD'], cwd=ROOT, text=True).strip()
    tracked = subprocess.check_output(['git', 'ls-files', '-z'], cwd=ROOT).decode('utf-8').split('\0')
    paths = set(p for p in tracked if p)
    paths.update(p.relative_to(ROOT).as_posix() for p in (ROOT / 'cloud').rglob('*')
                 if p.is_file() and not any(x in p.relative_to(ROOT / 'cloud').parts
                 for x in ('build', 'local_runs', '__pycache__')))
    paths.add('谷歌云端副本建立方案.md')
    excluded, included = [], []
    for name in sorted(paths):
        path = ROOT / name
        if name.startswith('tmp/') and not name.startswith('tmp/arxiv_2408.06515v2_source/'):
            excluded.append(dict(path=name, reason='temporary diagnostic output'))
        elif path.is_file():
            included.append(name)
    files = {name: dict(size_bytes=(ROOT/name).stat().st_size,
                       sha256=digest((ROOT/name).read_bytes())) for name in included}
    snapshot = revision[:7] + '-' + digest(json.dumps(files, sort_keys=True).encode())[:12]
    manifest = dict(schema=1, snapshot=snapshot, base_git_commit=revision,
        created_utc=datetime.now(timezone.utc).isoformat(),
        git_status=subprocess.check_output(['git', 'status', '--porcelain'], cwd=ROOT).decode('utf-8'),
        note='Current checkout, including the explicitly listed cloud/resume additions; not a clean export of base commit.',
        files=files, excluded_files=excluded,
        excluded_policies=['.git and editor/virtualenv caches', 'untracked files outside cloud and the supplied plan',
                           'cloud/build and cloud/local_runs', 'tmp except author source figures/text'],
        source_figure_exception='Keep original author PDF figures/text needed by the existing diagnostic scripts.')
    manifest_bytes = json.dumps(manifest, ensure_ascii=False, indent=2).encode('utf-8')
    archive = BUILD / f'pbh_{snapshot}.zip'
    with zipfile.ZipFile(archive, 'w', compression=zipfile.ZIP_DEFLATED) as bundle:
        for name in included:
            bundle.write(ROOT/name, 'project/'+name)
        bundle.writestr('project/SNAPSHOT_MANIFEST.json', manifest_bytes)
        bundle.writestr('project/cloud/source_changes.patch',
            subprocess.check_output(['git', 'diff', '--binary', 'HEAD'], cwd=ROOT))
    archive_hash = digest(archive.read_bytes())
    (BUILD/'manifest.json').write_bytes(manifest_bytes)
    with zipfile.ZipFile(archive) as bundle:
        if bundle.testzip() is not None:
            raise RuntimeError('ZIP CRC validation failed.')
        for name, record in files.items():
            if digest(bundle.read('project/'+name)) != record['sha256']:
                raise RuntimeError('Archive hash differs: '+name)
    markdown = f'''# PBH Colab：最新自适应积分与 cohort 方案 B

快照 `{snapshot}`，基于 Git `{revision[:7]}`，含可选 cohort 续算接口。
Google Drive 目录：`{DRIVE_NAME}`。选择 CPU 运行时，点击“运行全部”并授权 Drive。
科学计算在 `/content` 的隔离 Python 环境执行；检查点与结果保存回 Drive。

顺序：文件校验 → 安装锁定依赖 → 8 组固定初态跨平台比较 → B/full 与 B/GW-only。
默认每壳初始 128、每入晕边界 8 个样本，仅用于云端工作流验收，不能作论文定量结果。
保留严格联合抽样、Prada12-HMF、sesana_endpoint、200 Myr 环境步和 adaptive_log_a。
中断后再次“运行全部”可继续同一快照、同一配置下已保存的 cohort 状态。
不同配置会生成独立目录。`runtime_status.json` 记录真实运行状态。
'''
    setup = f'''from pathlib import Path
import datetime, hashlib, json, os, shutil, subprocess, sys, zipfile
from google.colab import drive
drive.mount('/content/drive')
DRIVE_ROOT = Path('/content/drive/MyDrive/{DRIVE_NAME}')
SNAPSHOT = {snapshot!r}
ARCHIVE_NAME = {archive.name!r}
ARCHIVE_SHA256 = {archive_hash!r}
RUN = DRIVE_ROOT / 'experiments' / SNAPSHOT
RUN.mkdir(parents=True, exist_ok=True)
def status(stage, **details):
    record = dict(stage=stage, snapshot=SNAPSHOT,
                  utc=datetime.datetime.now(datetime.timezone.utc).isoformat(), **details)
    temp = RUN / 'runtime_status.json.partial'
    temp.write_text(json.dumps(record, ensure_ascii=False, indent=2), encoding='utf-8')
    temp.replace(RUN / 'runtime_status.json')
    print(record)
status('drive_mounted', python=sys.version)
'''
    extract = '''status('verifying_snapshot')
try:
    source = DRIVE_ROOT / 'snapshots' / ARCHIVE_NAME
    work = Path('/content') / ('pbh_' + SNAPSHOT)
    work.mkdir(exist_ok=True)
    archive_local = work / ARCHIVE_NAME
    shutil.copyfile(source, archive_local)
    if hashlib.sha256(archive_local.read_bytes()).hexdigest() != ARCHIVE_SHA256:
        raise RuntimeError('快照校验失败，停止运行。')
    with zipfile.ZipFile(archive_local) as z:
        for member in z.infolist():
            if not (work / member.filename).resolve().is_relative_to(work.resolve()):
                raise RuntimeError('压缩包路径不合法。')
        z.extractall(work)
    PROJECT = work / 'project'
    manifest = json.loads((PROJECT / 'SNAPSHOT_MANIFEST.json').read_text(encoding='utf-8'))
    for name, record in manifest['files'].items():
        if hashlib.sha256((PROJECT / name).read_bytes()).hexdigest() != record['sha256']:
            raise RuntimeError('项目文件校验失败：' + name)
    if sys.version_info < (3, 12):
        raise RuntimeError('该版本 hmf 需要 Python >=3.12；请选择较新的 Colab 运行时。')
    status('snapshot_verified', files=len(manifest['files']))
except Exception as exc:
    status('failed', step='snapshot', error=str(exc))
    raise
'''
    install = '''def run_logged(command, log_name):
    local_log = work / log_name
    try:
        with local_log.open('w', encoding='utf-8') as log:
            process = subprocess.Popen([str(x) for x in command], cwd=PROJECT,
                stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True,
                env={**os.environ, 'OMP_NUM_THREADS':'1', 'OPENBLAS_NUM_THREADS':'1',
                     'MPLBACKEND':'Agg', 'PYTHONUNBUFFERED':'1'})
            for line in process.stdout:
                log.write(line)
                print(line, end='')
            code = process.wait()
        if code:
            raise RuntimeError(f'命令退出码 {code}，请查看 {log_name}')
    except BaseException as exc:
        status('failed_or_interrupted', log=log_name, error=str(exc))
        raise
    finally:
        if local_log.exists():
            shutil.copyfile(local_log, RUN / log_name)

status('creating_environment', launcher_revision='host_pip_v2')
ENV = work / 'venv'
run_logged([sys.executable, '-m', 'venv', '--without-pip', ENV], 'venv.log')
PYTHON = ENV / 'bin' / 'python'
run_logged([sys.executable, '-m', 'pip', '--version'], 'host_pip.log')
status('installing_dependencies', launcher_revision='host_pip_v2')
run_logged([sys.executable, '-m', 'pip', '--python', PYTHON,
            'install', '-r', PROJECT/'cloud/requirements-colab.txt'], 'install.log')
status('dependencies_installed')
'''
    comparison = '''status('running_cross_platform_comparison')
run_logged([PYTHON, PROJECT/'cloud/run_cloud.py', 'compare', '--output', RUN/'comparison'], 'comparison.log')
comparison = json.loads((RUN/'comparison/comparison.json').read_text(encoding='utf-8'))
if not comparison['passed']:
    raise RuntimeError('跨平台对照未通过，后续计算暂停。')
status('cross_platform_comparison_passed')
'''
    run = '''# 可编辑参数。初跑用于验证；增加样本前先检查运行时间与误差。
SEEDS = [20260923]
INITIAL_SAMPLES = 128
ACCRETION_SAMPLES = 8
GLOBAL_TIMESTEP_MYR = 200
ENTRY_MODEL = 'gw_aged'
if not json.loads((RUN/'comparison/comparison.json').read_text())['passed']:
    raise RuntimeError('必须先通过本次快照的跨平台对照。')
for seed in SEEDS:
    for mode in ('full', 'gw_only'):
        status('running_cohort', seed=seed, mode=mode)
        run_logged([PYTHON, PROJECT/'cloud/run_cloud.py', 'cohort', '--output', RUN/'jobs',
                    '--seed', seed, '--mode', mode, '--samples', INITIAL_SAMPLES,
                    '--accretion', ACCRETION_SAMPLES, '--timestep', GLOBAL_TIMESTEP_MYR,
                    '--entry', ENTRY_MODEL], f'cohort_{seed}_{mode}.log')
status('complete', note='Execution pilot completed; statistical convergence not claimed.')
'''
    display = '''from IPython.display import display, Image
for summary_path in sorted((RUN/'jobs').glob('*/summary.json')):
    summary = json.loads(summary_path.read_text(encoding='utf-8'))
    print({k:summary[k] for k in ('job','samples_drawn','inside_merged_weight',
                                 'outside_merged_weight','balance_max_abs')})
    display(Image(filename=str(summary_path.parent/'rates.png')))
print('结果已写回：', RUN)
'''
    cells = []
    for index, (kind, source) in enumerate([('markdown', markdown), ('code', setup),
        ('code', extract), ('code', install), ('code', comparison), ('code', run), ('code', display)]):
        cell = dict(cell_type=kind, id=f'pbh-cloud-{index}', metadata={}, source=source.splitlines(True))
        if kind == 'code':
            cell.update(execution_count=None, outputs=[])
        cells.append(cell)
    notebook = dict(cells=cells, metadata=dict(kernelspec=dict(display_name='Python 3',
        language='python', name='python3'), colab=dict(name='colab_bootstrap.ipynb'),
        language_info=dict(name='python')), nbformat=4, nbformat_minor=5)
    (BUILD/'colab_bootstrap.ipynb').write_text(json.dumps(notebook, ensure_ascii=False, indent=2), encoding='utf-8')
    receipt = dict(snapshot=snapshot, archive=archive.name, sha256=archive_hash,
                   md5=hashlib.md5(archive.read_bytes()).hexdigest(), size_bytes=archive.stat().st_size,
                   files=len(files), notebook='colab_bootstrap.ipynb')
    (BUILD/'bundle.json').write_text(json.dumps(receipt, indent=2), encoding='utf-8')
    print(json.dumps(receipt, indent=2))


if __name__ == '__main__':
    main()
