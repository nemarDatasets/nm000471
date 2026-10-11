"""src-literature-169 (OSF kpxw9 'Treder_et_al_2021', CC0-1.0; Treder et al. 2021 PNAS doi:10.1073/pnas.2114171118, "The hippocampus
as the switchboard between perception and memory") -> BIDS-iEEG DERIVATIVE (epoched iEEG).
Source: sharing/iEEG/s<NN>_{localiser,recall}_EEG_extraHIPP.mat (MAT v7.3, FieldTrip struct eeg: fsample, label 'channelNN', time{}, trial{}
(samples x ch, single), trialinfo{} structs). The *_recall_TFR_HIPP.mat files are time-frequency power (derived; not converted, in
sourcedata listing only) and sharing/EEG/ holds the scalp-EEG experiment (not iEEG; not converted).
Usage: python o5c_A1_convert_lit169.py <open500 record dir> <bids root>"""
import glob, json, os, re, sys
import numpy as np, h5py
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from o5c_A1_common import write_epochs, events_json, h5str, h5cellstr, wtsv, wjson, finalize_common, copy_code, copy_sourcedata, sha256_file

REC, OUT = sys.argv[1], sys.argv[2]
SRC = os.path.join(REC, "files", "sharing")
report = {"recordings": [], "skipped": [], "privacy_strings": []}


def tval(f, o):
    """trialinfo field -> scalar/str."""
    if isinstance(o, h5py.Group):
        return "struct"
    mc = o.attrs.get("MATLAB_class", b"")
    if mc == b"char":
        return h5str(f, o)
    a = np.asarray(o[()])
    if a.dtype == object or mc == b"cell":
        refs = a.ravel()
        return tval(f, f[refs[0]]) if refs.size == 1 else ";".join(str(tval(f, f[r])) for r in refs)
    if a.dtype.kind in "ui" and a.size and "MATLAB_empty" in o.attrs:
        return "n/a"
    if a.size == 1:
        v = float(a.ravel()[0])
        if not np.isfinite(v):
            return "n/a" if np.isnan(v) else ("Inf" if v > 0 else "-Inf")
        return int(v) if v == int(v) else v
    return "n/a" if a.size == 0 else ";".join(str(x) for x in a.ravel().tolist()[:20])


scans = {}
files = sorted(glob.glob(os.path.join(SRC, "iEEG", "s*_localiser_EEG_extraHIPP.mat")) + glob.glob(os.path.join(SRC, "iEEG", "s*_recall_EEG_extraHIPP.mat")))
for p in files:
    m = re.match(r"s(\d+)_(localiser|recall)_EEG_extraHIPP\.mat$", os.path.basename(p))
    sid, task = int(m[1]), m[2]
    sub = f"sub-{sid:02d}"; D = os.path.join(OUT, sub, "ieeg"); stem = f"{sub}_task-{task}"
    with h5py.File(p, "r") as f:
        e = f["eeg"]
        fs = float(np.asarray(e["fsample"][()]).ravel()[0])
        labels = h5cellstr(f, e["label"])
        T = np.asarray(e["trial"][()]).ravel(); TM = np.asarray(e["time"][()]).ravel(); TI = np.asarray(e["trialinfo"][()]).ravel()
        assert T.size == TM.size == TI.size, (T.size, TM.size, TI.size)
        trials, times, rows = [], [], []
        keys = None
        for k in range(T.size):
            x = np.asarray(f[T[k]][()])  # h5 (samples, ch)
            if x.shape[1] != len(labels) and x.shape[0] == len(labels):
                x = x.T
            t = np.asarray(f[TM[k]][()]).ravel()
            assert x.shape[0] == t.size, (x.shape, t.shape)
            g = f[TI[k]]
            if keys is None:
                keys = list(g.keys())
            rows.append([tval(f, g[kk]) if kk in g else "n/a" for kk in keys])
            trials.append(x.astype(np.float64)); times.append(t)
        for s_ in h5cellstr(f, e["label"]):
            pass
    fs_r = round(fs, 3)
    r = write_epochs(D, stem, trials, fs_r, labels, ["n/a"] * len(labels), ["SEEG"] * len(labels), {
        "TaskName": task,
        "TaskDescription": ("Localiser: object or scene cue words shown while the patient performed a perceptual categorisation (article); trialinfo: StimType (1 object cue, 2 scene cue per the authors' script PNAS_F3_B.m), Response, RT, Accuracy, RunNumber, TrialNumber." if task == "localiser" else
                            "Associative recall: after learning word-picture (object/scene) pairs, cue words were shown and the patient recalled the associated picture (article); trialinfo: Memory (R remembered / F forgot), ObjectProbe, Response, RT, ExpPhase, RunNumber, TrialNumber."),
        "Instructions": "n/a (see article)",
        "PowerLineFrequency": 50,
        "SoftwareFilters": {"AuthorsPreprocessing": {"Description": "Released as FieldTrip epochs at 200 Hz (single precision); preprocessing as described in the article (not in the release). Values written unchanged."}},
        "HardwareFilters": "n/a",
        "iEEGReference": "n/a (not stated in the release)",
        "iEEGPlacementScheme": "intracranial depth electrodes implanted for epilepsy monitoring; the released channels are labelled channel01..N (anonymised; file name 'extraHIPP'). No coordinates or anatomical labels are released for these channels.",
        "ElectricalStimulation": False, "Manufacturer": "n/a",
    }, keys, rows, ch_extra=[{"description": "released label %s (anonymised contact name; location not released)" % l} for l in labels],
        comment=f"o5c_A1_convert_lit169.py: {os.path.basename(p)} eeg.trial back-to-back, float32", trial_times=times)
    events_json(os.path.join(D, stem + "_events.json"), {k: {"Description": f"eeg.trialinfo.{k} as released"} for k in keys})
    r.update(sub=sub, task=task, file=os.path.basename(p), fs_raw=fs, trialinfo_keys=keys)
    report["recordings"].append(r)
    scans.setdefault(sub, []).append([f"ieeg/{stem}_ieeg.vhdr", "n/a"])
    print(stem, json.dumps({k: v for k, v in r.items() if k != "trialinfo_keys"}), flush=True)
for sub, rows in scans.items():
    wtsv(os.path.join(OUT, sub, f"{sub}_scans.tsv"), ["filename", "acq_time"], sorted(rows))
subs = sorted(scans)
wtsv(os.path.join(OUT, "participants.tsv"), ["participant_id", "age", "sex", "source_id"], [[s, "n/a", "n/a", "s%02d" % int(s[4:])] for s in subs])
wjson(os.path.join(OUT, "participants.json"), {"age": {"Description": "not in the release"}, "sex": {"Description": "not in the release"}, "source_id": {"Description": "subject code in the release file names"}})
from o5c_A1_common import add_electrodes_na
add_electrodes_na(OUT, "No electrode coordinates or anatomical labels are released for the iEEG channels (anonymised names channel01..N; the file names say 'extraHIPP'). x, y, z = n/a.", spatial_ref="n/a (no coordinates released)")
# sourcedata: the iEEG epoch files + scripts (the TFR power files and the scalp-EEG folder are listed, not copied)
sd = []
for p in sorted(glob.glob(os.path.join(SRC, "iEEG", "*_EEG_extraHIPP.mat")) + glob.glob(os.path.join(SRC, "*.m")) + glob.glob(os.path.join(SRC, "auxiliary", "*"))):
    sd.append((p, os.path.relpath(p, SRC)))
SD = os.path.join(OUT, "sourcedata", "osf-kpxw9")
copy_sourcedata(sd, SD, scrub=False)
listed = [[os.path.relpath(p, SRC), os.path.getsize(p)] for p in sorted(glob.glob(os.path.join(SRC, "iEEG", "*_TFR_HIPP.mat")) + glob.glob(os.path.join(SRC, "EEG", "*")))]
wtsv(os.path.join(SD, "NOT_COPIED_release_files.tsv"), ["path", "bytes"], listed)
copy_code(OUT, [__file__])
Rr = report["recordings"]; hrs = sum(r["duration_s"] for r in Rr) / 3600
finalize_common(OUT,
    name="Treder2021 hippocampal switchboard: epoched intracranial EEG during object/scene localiser and associative recall (extra-hippocampal channels)",
    dtype="derivative", license_="CC0-1.0",
    authors=["Matthias S. Treder", "Ian Charest", "Sebastian Michelmann", "María Carmen Martín-Buro", "Frédéric Roux", "Fernando Carceller-Benito", "Arturo Ugalde-Canitrot", "David T. Rollings", "Vijay Sawlani", "Ramesh Chelvarajah", "Maria Wimber", "Simon Hanslmayr", "Bernhard P. Staresina"],
    refs=["https://doi.org/10.1073/pnas.2114171118", "https://doi.org/10.1101/2020.05.20.104539", "https://osf.io/kpxw9/"],
    howto="Cite Treder MS, Charest I, Michelmann S, et al. (2021) The hippocampus as the switchboard between perception and memory. PNAS 118(50):e2114171118, doi:10.1073/pnas.2114171118, and the OSF release https://osf.io/kpxw9/.",
    keywords=["intracranial EEG", "hippocampus", "episodic memory", "associative recall", "object", "scene", "epilepsy"],
    source=[{"URL": "https://osf.io/kpxw9/", "Version": "OSF node kpxw9 'Treder_et_al_2021' (files as downloaded 2026-10-11; sha256 in sourcedata/osf-kpxw9/DEIDENTIFICATION_MANIFEST.tsv)"}],
    generated_by=[{"Name": "o5c_A1_convert_lit169.py", "Description": "Repackaging of the FieldTrip epoch structs (eeg.trial) into BrainVision float32, trials back-to-back, with eeg.trialinfo as events. Code in code/."}],
    ethics=["Collected under the source study's approval (see the article, doi:10.1073/pnas.2114171118)."],
    readme=f"""# Treder et al. 2021, "The hippocampus as the switchboard between perception and memory": intracranial EEG epochs

Epoched intracranial EEG from {len(subs)} patients ({len(Rr)} recordings, {hrs:.2f} h of epochs, 200 Hz) during an object/scene
localiser task and an associative recall task, released by the authors on OSF (https://osf.io/kpxw9/, CC0) with the MATLAB
scripts that reproduce the PNAS figures.

## What is converted (DatasetType derivative)
- `sub-XX/ieeg/sub-XX_task-localiser_ieeg.vhdr` and `..._task-recall_ieeg.vhdr`: the released `s<NN>_<task>_EEG_extraHIPP.mat`
  FieldTrip epochs (single precision, 200 Hz), written back-to-back with one `New Segment` marker per trial; `events.tsv` has
  each trial's position in the file, the authors' time axis start (-2 s) and every `trialinfo` field as released.
- Channel names are the released anonymised labels `channel01..N`; no coordinates or anatomical labels are released for them
  (`electrodes.tsv` with n/a coordinates). Per the file names these are the extra-hippocampal channels; the hippocampal data of
  the release are time-frequency power only (`*_recall_TFR_HIPP.mat`, derived, not converted; listed in
  `sourcedata/osf-kpxw9/NOT_COPIED_release_files.tsv` together with the scalp-EEG experiment files of `sharing/EEG/`).
- Signal units are not stated in the release (`units` = n/a; values unchanged). Preprocessing is described in the article.
- `sourcedata/osf-kpxw9/`: the converted MAT files, the authors' scripts and auxiliary files, unchanged.
- Age and sex are not in the release.

## Privacy
Channel names are anonymised by the authors; files hold trial numbers, responses and reaction times only.
""",
    changes_note="BIDS conversion of the OSF kpxw9 release (Treder et al. 2021) for NEMAR (o5c lane A1).")
json.dump(report, open(os.path.join(OUT, "code", "conversion_report.json"), "w"), indent=1, default=str)
print(json.dumps({"recordings": len(Rr), "subjects": len(subs), "hours": hrs, "rt_all": all(r["roundtrip_ok"] for r in Rr), "channels": [r["n_channels"] for r in Rr]}, indent=1))
print("CONVERT_DONE")
