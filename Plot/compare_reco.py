#!/usr/bin/env python3
import argparse
import os
import numpy as np

def main():
    try:
        import uproot
    except Exception:
        print('Missing dependency: uproot. Install with `pip install uproot awkward numpy matplotlib`')
        raise

    parser = argparse.ArgumentParser(description='Compare reconstruction methods')
    parser.add_argument('rootfile', help='Input ROOT file')
    parser.add_argument('--region', default='ttll_fully_matched', help='Region suffix used in branch names')
    parser.add_argument('--outdir', default='analysis_plots', help='Output directory')
    parser.add_argument('--nbins', type=int, default=100, help='Number of bins for mass histograms')
    parser.add_argument('--xmax', type=float, default=2000.0, help='Max m(ttbar) x-axis in GeV (will cap hist range)')
    args = parser.parse_args()

    path = args.rootfile
    outdir = args.outdir
    os.makedirs(outdir, exist_ok=True)

    f = uproot.open(path)
    if 'NOSYS' not in f:
        keys = [k for k in f.keys()]
        if not keys:
            raise RuntimeError('No keys found in ROOT file')
        tree = f[keys[0]]
    else:
        tree = f['NOSYS']

    ks = [k if isinstance(k, str) else k.decode('utf-8') for k in tree.keys()]
    region = args.region

    # try to prefer TTree named 'reco' or a directory named after the region
    top_keys = [k if isinstance(k, str) else k.decode('utf-8') for k in f.keys()]
    if 'reco' in top_keys:
        tree = f['reco']
    elif region in top_keys:
        tree = f[region]
    # recompute ks for the selected tree
    ks = [k if isinstance(k, str) else k.decode('utf-8') for k in tree.keys()]

    # helper: rebin histogram counts to target edges by overlap
    def rebin_to_target(counts, edges, target_edges):
        new_counts = np.zeros(len(target_edges) - 1, dtype=float)
        for i in range(len(counts)):
            a0 = edges[i]
            a1 = edges[i+1]
            width = a1 - a0
            if width <= 0:
                continue
            for j in range(len(new_counts)):
                b0 = target_edges[j]
                b1 = target_edges[j+1]
                overlap = max(0.0, min(a1, b1) - max(a0, b0))
                if overlap > 0:
                    new_counts[j] += counts[i] * (overlap / width)
        return new_counts

    # function to find an existing suffix for a requested multiplicity
    def find_suffix_for_mult(mult):
        cand_list = []
        if isinstance(mult, str) and mult.endswith('+'):
            nmin = mult[:-1]
            cand_list = [
                f'{nmin}jets', f'{nmin}+jets', f'{nmin}plus', f'{nmin}plusjets',
                f'ge{nmin}jets', f'>={nmin}jets'
            ]
        else:
            n = str(mult)
            cand_list = [f'{n}jets', f'{n}jet', f'{n}_jets', f'{n}_jet']
        for cand in cand_list:
            for k in ks:
                # match case-insensitive ending like _{region}_{cand}
                if k.lower().endswith(f'_{region}_{cand}'.lower()):
                    return cand
        return None

    # process a single jet multiplicity (integer or '5+') and handle both
    # histogram outputs and event-level ntuple branches
    def process_multiplicity(mult):
        import math
        import awkward as ak

        # helper: build candidate branch names for correctness per method
        # Prefer the per-method boolean correctness branches (PT/CHI2/DELTAR/HYPER),
        # falling back to legacy pairing scalar branches if needed.
        method_patterns = {
            'baseline': [
                f'PT_correct_{region}_{mult}', f'PT_correct_{region}', 'PT_correct_NOSYS', 'PT_correct',
                f'baseline_pairing_correct_{region}_{mult}', f'baseline_pairing_correct_{region}', 'baseline_pairing_correct', 'pairing_correct_NOSYS', 'pairing_correct'
            ],
            'chisquared': [
                f'CHI2_correct_{region}_{mult}', f'CHI2_correct_{region}', 'CHI2_correct_NOSYS', 'CHI2_correct',
                f'chisquared_pairing_correct_all_{region}_{mult}', f'chisquared_pairing_correct_all_{region}', 'chisquared_pairing_correct_all', 'pairing_chi2_correct_all_NOSYS', 'pairing_chi2_correct_all'
            ],
            'deltaR': [
                f'DELTAR_correct_{region}_{mult}', f'DELTAR_correct_{region}', 'DELTAR_correct_NOSYS', 'DELTAR_correct',
                f'deltaR_pairing_correct_{region}_{mult}', f'deltaR_pairing_correct_{region}', 'deltaR_pairing_correct', 'pairing_deltaR_correct_NOSYS', 'pairing_deltaR_correct'
            ],
            'hyper': [
                f'HYPER_correct_{region}_{mult}', f'HYPER_correct_{region}', 'HYPER_correct_NOSYS', 'HYPER_correct'
            ]
        }

        # collect correctness info — histograms (TH1) or event-level awkward arrays
        found = {}
        for m, pats in method_patterns.items():
            candidate = None
            for p in pats:
                if not p:
                    continue
                if p in ks:
                    candidate = p
                    break
                if (p + ';1') in ks:
                    candidate = p + ';1'
                    break
            if candidate is None:
                continue
            obj = tree[candidate]
            # TH1 objects provide to_numpy; TTree branches return arrays
            try:
                if hasattr(obj, 'to_numpy'):
                    counts, edges = obj.to_numpy()
                    found[m] = ('hist', counts, edges)
                else:
                    # use awkward arrays for event-level branches
                    arr = obj.array(library='ak')
                    found[m] = ('ak', arr)
            except Exception:
                try:
                    # fallback to numpy arrays
                    arr = obj.array(library='np')
                    found[m] = ('array', np.asarray(arr))
                except Exception:
                    continue

        if not found:
            print('No method correctness branches found for region', region, 'multiplicity', mult)
            return

        eff = {}
        err = {}
        totals = {}
        # If we have event-level arrays, we need to filter by multiplicity
        def get_njets_array():
            # Prefer multiplicity inferred from the size of a jet_pt-like jagged branch
            # (use ak.num or per-event .size()/len), falling back to explicit n_jets
            jet_cands = ['jet_pt_NOSYS', 'jet_pt', 'Jet_pt', 'Jet_pt_NOSYS', 'jet_TLV_NOSYS', 'jet_TLV', 'jets_pt', 'jets_pt_NOSYS']
            for cand in jet_cands:
                if cand in ks or (cand + ';1') in ks:
                    k = cand if cand in ks else cand + ';1'
                    # Try awkward jagged arrays first and use ak.num (equivalent to .size())
                    try:
                        arr_ak = tree[k].array(library='ak')
                        try:
                            nums = ak.num(arr_ak)
                            return np.asarray(ak.to_numpy(nums), dtype=int)
                        except Exception:
                            # Last-resort: iterate and compute length/size per event
                            counts = []
                            for elem in arr_ak:
                                try:
                                    counts.append(len(elem))
                                except Exception:
                                    try:
                                        counts.append(int(getattr(elem, 'size', 0) or 0))
                                    except Exception:
                                        counts.append(0)
                            return np.asarray(counts, dtype=int)
                    except Exception:
                        # Fallback to numpy object arrays where each entry is an iterable
                        try:
                            arr_np = tree[k].array(library='np')
                            if isinstance(arr_np, np.ndarray) and arr_np.dtype == object:
                                return np.array([len(x) if x is not None else 0 for x in arr_np], dtype=int)
                            # If we get here, arr_np may already be a per-event integer array; return as-is
                            return np.asarray(arr_np, dtype=int)
                        except Exception:
                            continue

            # Fall back to explicit n_jets-like branches if no jet_pt-like branch found
            nj_cands = ['n_jets_NOSYS', 'n_jets', 'nJets', 'nJets_good_NOSYS', 'nJets_good', 'n_jets_good_NOSYS', 'n_jets_good']
            for cand in nj_cands:
                if cand in ks or (cand + ';1') in ks:
                    k = cand if cand in ks else cand + ';1'
                    try:
                        arr_np = tree[k].array(library='np')
                        return np.asarray(arr_np, dtype=int)
                    except Exception:
                        try:
                            arr_ak = tree[k].array(library='ak')
                            try:
                                return np.asarray(ak.to_numpy(arr_ak), dtype=int)
                            except Exception:
                                return np.asarray([int(x) if x is not None else 0 for x in arr_ak], dtype=int)
                        except Exception:
                            continue
            return None

        njets_arr = get_njets_array()

        def multiplicity_mask(njets, mult_label):
            if njets is None:
                return None
            arr = np.asarray(njets)
            s = str(mult_label)
            if s.endswith('+'):
                try:
                    nmin = int(s[:-1])
                    return arr >= nmin
                except Exception:
                    return np.ones(arr.shape[0], dtype=bool)
            try:
                n = int(s.lstrip('_').split('j')[0])
                return arr == n
            except Exception:
                return np.ones(arr.shape[0], dtype=bool)

        for m, info in found.items():
            if info[0] == 'hist':
                counts, edges = info[1], info[2]
                centers = 0.5 * (edges[:-1] + edges[1:])
                idx_true = int(np.argmin(np.abs(centers - 1.0)))
                n_true = float(counts[idx_true])
                total = float(np.sum(counts))
            else:
                # event-level arrays
                if info[0] == 'ak':
                    arr_ak = info[1]
                    # convert full array to numpy first (robust for scalars)
                    try:
                        vals_all = ak.to_numpy(arr_ak)
                    except Exception:
                        # best-effort conversion
                        vals_all = np.asarray([x if x is not None else 0 for x in arr_ak])

                    # apply multiplicity mask if available (njets_arr is returned as numpy by get_njets_array)
                    try:
                        if njets_arr is not None:
                            mask = multiplicity_mask(njets_arr, mult)
                            vals = np.asarray(vals_all)[mask]
                        else:
                            vals = np.asarray(vals_all)

                        vals = np.asarray(vals)
                        # determine truth count: booleans or numeric flags where 1 means correct
                        #vals = np.asarray(vals)

                        # keep only valid events (0 or 1)
                        valid_mask = (vals == 0) | (vals == 1)
                        vals = vals[valid_mask]

                        total = float(vals.size)

                        if total > 0:
                            n_true = float(np.count_nonzero(vals == 1))
                        else:
                            n_true = 0.0
                    except Exception:
                        n_true = 0.0
                        total = 0.0
                else:
                    arr = info[1]
                    vals = np.asarray(arr).astype(bool)
                    # if multiplicity available, need to filter — but numpy object arrays are hard; skip filtering
                    n_true = float(np.count_nonzero(vals))
                    total = float(len(vals))
            e = n_true / total if total > 0 else 0.0
            ee = math.sqrt(e * (1.0 - e) / total) if total > 0 else 0.0
            eff[m] = e
            err[m] = ee
            totals[m] = int(total)

        print('Region:', region, 'multiplicity:', mult)
        for m in eff:
            print(f'  {m}: {eff[m]:.4f} ± {err[m]:.4f} (n_true~={int(eff[m]*totals[m])}, N={totals[m]})')

        # save bar plot for efficiencies for this multiplicity
        try:
            import matplotlib
            matplotlib.use('Agg')
            import matplotlib.pyplot as plt
        except Exception:
            print('matplotlib not available; install with pip to save plots')
            return

        methods_sorted = list(eff.keys())
        y = [eff[m] for m in methods_sorted]
        yerr = [err[m] for m in methods_sorted]
        fig, ax = plt.subplots(figsize=(6,4))
        palette = ['#1f77b4','#ff7f0e','#2ca02c','#d62728']
        ax.bar(methods_sorted, y, yerr=yerr, capsize=5, color=palette[:len(methods_sorted)])
        ax.set_ylim(0,1)
        ax.set_ylabel('Correct reconstruction efficiency')
        suffix_str = str(mult)
        suffix_safe = suffix_str.replace('+','plus')
        ax.set_title(f'Efficiency per method — {region} ({suffix_str})')
        for i, v in enumerate(y):
            ax.text(i, v + 0.02, f'{v:.3f}', ha='center')
        fout = os.path.join(outdir, f'efficiency_{region}_{suffix_safe}.png')
        fig.savefig(fout, dpi=150)
        print('Saved efficiency plot to', fout)
        plt.close(fig)

        # m(ttbar) distributions: plot both NUFLOWS and NUWEIGHT variants.
        method_vars = [
            'mttbar_HYPER_NUFLOWS', 'mttbar_CHI2_NUFLOWS', 'mttbar_PT_NUFLOWS', 'mttbar_DELTAR_NUFLOWS',
            'mttbar_HYPER_NUWEIGHT', 'mttbar_CHI2_NUWEIGHT', 'mttbar_PT_NUWEIGHT', 'mttbar_DELTAR_NUWEIGHT'
        ]
        legacy_nuflows_map = {
            'mttbar_HYPER_NUFLOWS': 'mttbar_HYPER',
            'mttbar_CHI2_NUFLOWS': 'mttbar_CHI2',
            'mttbar_PT_NUFLOWS': 'mttbar_PT',
            'mttbar_DELTAR_NUFLOWS': 'mttbar_DELTAR',
        }
        # collect per-method histograms/arrays
        method_data = {}
        for mv in method_vars:
            # look for multiple candidate names (with region/suffix, NOSYS, or bare)
            pats = [f'{mv}_{region}_{mult}', f'{mv}_{region}', f'{mv}_NOSYS', mv]
            legacy = legacy_nuflows_map.get(mv)
            if legacy is not None:
                pats.extend([f'{legacy}_{region}_{mult}', f'{legacy}_{region}', f'{legacy}_NOSYS', legacy])
            cand = None
            for p in pats:
                if p in ks:
                    cand = p
                    break
                if (p + ';1') in ks:
                    cand = p + ';1'
                    break
            if not cand:
                continue
            obj = tree[cand]
            try:
                if hasattr(obj, 'to_numpy'):
                    counts, edges = obj.to_numpy()
                    method_data[mv] = ('hist', counts, edges)
                else:
                    arr = obj.array(library='ak')
                    method_data[mv] = ('ak', arr)
            except Exception:
                try:
                    arr = obj.array(library='np')
                    method_data[mv] = ('array', np.asarray(arr))
                except Exception:
                    continue

        # truth per-multiplicity
        # truth variable: prefer mttbar_truth if available, otherwise llbb_mass / invariant_mass_llbb
        truth_data = None
        truth_pats = [
            f'mttbar_TRUTH_{region}_{mult}', f'mttbar_TRUTH_{region}', 'mttbar_TRUTH_NOSYS', 'mttbar_TRUTH',
            f'mttbar_truth_{region}_{mult}', f'mttbar_truth_{region}', 'mttbar_truth_NOSYS', 'mttbar_truth'
        ]
        for tpat in truth_pats:
            cand = None
            if tpat in ks:
                cand = tpat
            elif (tpat + ';1') in ks:
                cand = tpat + ';1'
            if not cand:
                continue
            obj = tree[cand]
            try:
                if hasattr(obj, 'to_numpy'):
                    truth_counts, truth_edges = obj.to_numpy()
                    truth_data = ('hist', truth_counts, truth_edges)
                else:
                    arr = obj.array(library='ak')
                    truth_data = ('ak', arr)
                break
            except Exception:
                try:
                    truth_arr = np.asarray(obj.array(library='np'))
                    truth_data = ('array', truth_arr)
                    break
                except Exception:
                    continue

        if not method_data and truth_data is None:
            print('No per-method mttbar or truth branches available for plotting for region', region, 'multiplicity', mult)
            return

        # determine unified bin edges across methods and truth
        lows = []
        highs = []
        for info in method_data.values():
            if info[0] == 'hist':
                _, edges = info[1], info[2]
                lows.append(edges[0]); highs.append(edges[-1])
            elif info[0] == 'ak':
                arr = info[1]
                try:
                    npv = ak.to_numpy(arr)
                    if npv.size > 0:
                        lows.append(float(np.nanmin(npv))); highs.append(float(np.nanmax(npv)))
                except Exception:
                    continue
            else:
                arr = info[1]
                try:
                    if arr.size > 0:
                        lows.append(float(np.nanmin(arr))); highs.append(float(np.nanmax(arr)))
                except Exception:
                    continue
        if truth_data is not None:
            if truth_data[0] == 'hist':
                lows.append(truth_data[2][0]); highs.append(truth_data[2][-1])
            elif truth_data[0] == 'ak':
                try:
                    npv = ak.to_numpy(truth_data[1])
                    if npv.size > 0:
                        lows.append(float(np.nanmin(npv))); highs.append(float(np.nanmax(npv)))
                except Exception:
                    pass
            else:
                arr = truth_data[1]
                try:
                    if arr.size > 0:
                        lows.append(float(np.nanmin(arr))); highs.append(float(np.nanmax(arr)))
                except Exception:
                    pass
        lo = min(lows) if lows else 0.0
        hi = max(highs) if highs else 1.0
        # cap hi at user-requested maximum (in GeV)
        try:
            hi = min(hi, float(args.xmax))
        except Exception:
            pass
        if lo == hi:
            lo = 0.0
            hi = max(1.0, hi)
        target_edges = np.linspace(lo, hi, args.nbins + 1)

        # build rebinned counts for each method
        method_counts = {}
        for mv, info in method_data.items():
            if info[0] == 'hist':
                counts, edges = info[1], info[2]
                method_counts[mv] = rebin_to_target(counts, edges, target_edges)
            elif info[0] == 'ak':
                try:
                    arr = ak.to_numpy(info[1])
                    if njets_arr is not None:
                        mask = multiplicity_mask(njets_arr, mult)
                        arr = arr[mask]
                    cts, _ = np.histogram(arr, bins=target_edges)
                    method_counts[mv] = cts
                except Exception:
                    method_counts[mv] = np.zeros(len(target_edges) - 1)
            else:
                try:
                    arr = np.asarray(info[1])
                    cts, _ = np.histogram(arr, bins=target_edges)
                    method_counts[mv] = cts
                except Exception:
                    method_counts[mv] = np.zeros(len(target_edges) - 1)

        if truth_data is not None:
            if truth_data[0] == 'hist':
                truth_counts = rebin_to_target(truth_data[1], truth_data[2], target_edges)
            elif truth_data[0] == 'ak':
                try:
                    tarr = ak.to_numpy(truth_data[1])
                    if njets_arr is not None:
                        mask = multiplicity_mask(njets_arr, mult)
                        tarr = tarr[mask]
                    truth_counts, _ = np.histogram(tarr, bins=target_edges)
                except Exception:
                    truth_counts = np.zeros(len(target_edges) - 1)
            else:
                try:
                    truth_counts, _ = np.histogram(truth_data[1], bins=target_edges)
                except Exception:
                    truth_counts = np.zeros(len(target_edges) - 1)
        else:
            truth_counts = np.zeros(len(target_edges) - 1)

        bin_centers = 0.5 * (target_edges[:-1] + target_edges[1:])
        # scale masses for plotting (divide by 1e3 so axes are in TeV)
        scale_mass = 1000.0
        bin_centers = bin_centers / scale_mass
        lo_plot = lo / scale_mass
        hi_plot = hi / scale_mass
        # normalized densities
        method_density = {mv: (cts / np.sum(cts) if np.sum(cts) > 0 else cts) for mv, cts in method_counts.items()}
        truth_density = truth_counts / np.sum(truth_counts) if np.sum(truth_counts) > 0 else truth_counts

        import matplotlib.pyplot as plt

        colors = ['#1f77b4','#ff7f0e','#2ca02c','#d62728']
        style_map = {'baseline': '-', 'chisquared': '--', 'deltaR': '-.', 'hyper': ':'}

        def mv_to_key(mv):
            l = mv.lower()
            if 'chi' in l:
                return 'chisquared'
            if 'delt' in l or 'delta' in l:
                return 'deltaR'
            if 'pt' in l:
                return 'baseline'
            if 'hyper' in l:
                return 'hyper'
            return None

        # Compose efficiency text box
        eff_lines = []
        for k in eff:
            eff_lines.append(f"{k}: {eff[k]:.3f} ± {err[k]:.3f}")
        eff_text = '\n'.join(eff_lines)
        bbox = dict(boxstyle='round', facecolor='white', alpha=0.8, edgecolor='none')

        # Counts plot with ratio subplot
        fig2, (ax2, ax2r) = plt.subplots(2, 1, sharex=True, figsize=(8, 8), gridspec_kw={'height_ratios': [3, 1]})
        fig2.subplots_adjust(hspace=0.05)
        for i, (mv, cts) in enumerate(method_counts.items()):
            key = mv_to_key(mv)
            label = mv.replace('mttbar_', '')
            if key in eff:
                label = f"{label} (eff={eff[key]:.3f})"
            ls = style_map.get(key, '-')
            ax2.step(bin_centers, cts, where='mid', label=label, color=colors[i % len(colors)], linestyle=ls)
        if truth_counts.sum() > 0:
            ax2.fill_between(bin_centers, truth_counts, step='mid', color='0.15', alpha=0.35, label='Truth (llbb_mass)')
            ax2.step(bin_centers, truth_counts, where='mid', label='Truth (llbb_mass)', color='k', linewidth=1.0)
        ax2.set_ylabel('Counts')
        ax2.set_title(f'm(ttbar) per-method vs truth — {region} ({suffix_str}) (counts)')
        ax2.legend()
        ax2.text(0.98, 0.95, eff_text, transform=ax2.transAxes, fontsize=8, va='top', ha='right', bbox=bbox)

        # Ratio to truth (counts)
        if truth_counts.sum() > 0:
            for i, (mv, cts) in enumerate(method_counts.items()):
                ratio = np.full_like(truth_counts, np.nan)
                mask = truth_counts > 0
                if np.any(mask):
                    ratio[mask] = cts[mask] / truth_counts[mask]
                ax2r.step(bin_centers, ratio, where='mid', color=colors[i % len(colors)])
            ax2r.axhspan(0.8, 1.2, color='gray', alpha=0.12)
            ax2r.set_ylim(0.5, 1.5)
        ax2r.axhline(1.0, color='k', linestyle='--', linewidth=1)
        ax2r.set_ylabel('Ratio')
        ax2r.set_xlabel('Mass [TeV]')
        ax2.grid(True, linestyle=':', linewidth=0.5, alpha=0.6)

        fout_norm_ratio = os.path.join(outdir, f'mttbar_methods_vs_truth_{region}_{suffix_safe}_norm_ratio.png')
        fig.savefig(fout_norm_ratio, dpi=150, bbox_inches='tight')
        print('Saved normalized per-method overlay with ratio to', fout_norm_ratio)
        plt.close(fig)

        # Normalized plot with ratio subplot
        fig, (ax, axr) = plt.subplots(2, 1, sharex=True, figsize=(8, 8), gridspec_kw={'height_ratios': [3, 1]})
        fig.subplots_adjust(hspace=0.05)
        for i, (mv, dens) in enumerate(method_density.items()):
            key = mv_to_key(mv)
            label = mv.replace('mttbar_', '')
            if key in eff:
                label = f"{label} (eff={eff[key]:.3f})"
            ls = style_map.get(key, '-')
            ax.step(bin_centers, dens, where='mid', label=label, color=colors[i % len(colors)], linestyle=ls)
        if truth_counts.sum() > 0:
            ax.fill_between(bin_centers, truth_density, step='mid', color='0.15', alpha=0.4, label='Truth (llbb_mass)')
            ax.step(bin_centers, truth_density, where='mid', color='k', linewidth=1.0)
        ax.set_ylabel('Normalized density')
        ax.set_title(f'm(ttbar) per-method vs truth — {region} ({suffix_str}) (normalized)')
        ax.legend()
        ax.text(0.98, 0.95, eff_text, transform=ax.transAxes, fontsize=8, va='top', ha='right', bbox=bbox)

        # Ratio to truth (normalized)
        if truth_density.sum() > 0:
            for i, (mv, dens) in enumerate(method_density.items()):
                ratio = np.full_like(truth_density, np.nan)
                mask = truth_density > 0
                if np.any(mask):
                    ratio[mask] = dens[mask] / truth_density[mask]
                axr.step(bin_centers, ratio, where='mid', color=colors[i % len(colors)])
            # highlight acceptable band and set reasonable y-limits for ratios
            axr.axhspan(0.8, 1.2, color='gray', alpha=0.12)
            axr.set_ylim(0.5, 1.5)
        axr.axhline(1.0, color='k', linestyle='--', linewidth=1)
        axr.set_ylabel('Ratio')
        axr.set_xlabel('Mass [TeV]')
        ax.grid(True, linestyle=':', linewidth=0.5, alpha=0.6)

        fout_counts_ratio = os.path.join(outdir, f'mttbar_methods_vs_truth_{region}_{suffix_safe}_counts_ratio.png')
        fig2.savefig(fout_counts_ratio, dpi=150, bbox_inches='tight')
        print('Saved counts per-method overlay with ratio to', fout_counts_ratio)
        plt.close(fig2)

        # Add a zoomed x-axis view for m(ttbar) (0..2000 GeV) if the data extend beyond that
        try:
            zoom_hi_limit = float(args.xmax)
            zoom_lo_plot = max(lo, 0.0) / scale_mass
            zoom_hi_plot = min(hi, zoom_hi_limit) / scale_mass
            if zoom_lo_plot < zoom_hi_plot:
                mask = (bin_centers >= zoom_lo_plot) & (bin_centers <= zoom_hi_plot)
                if np.any(mask):
                    zoom_centers = bin_centers[mask]
                    figz, (axz, axzr) = plt.subplots(2, 1, sharex=True, figsize=(8, 8), gridspec_kw={'height_ratios': [3, 1]})
                    figz.subplots_adjust(hspace=0.05)
                    for i, (mv, dens) in enumerate(method_density.items()):
                        dens_sub = dens[mask]
                        label = mv.replace('mttbar_', '')
                        key = mv_to_key(mv)
                        if key in eff:
                            label = f"{label} (eff={eff[key]:.3f})"
                        axz.step(zoom_centers, dens_sub, where='mid', label=label, color=colors[i % len(colors)])
                    if truth_counts.sum() > 0:
                        truth_sub = truth_density[mask]
                        axz.fill_between(zoom_centers, truth_sub, step='mid', color='0.15', alpha=0.35, label='Truth (llbb_mass)')
                        axz.step(zoom_centers, truth_sub, where='mid', label='Truth', color='k', linewidth=1.0)
                    axz.set_ylabel('Normalized density')
                    axz.set_title(f'm(ttbar) per-method vs truth — {region} ({suffix_str}) (zoom {zoom_lo_plot:.2f}-{zoom_hi_plot:.2f} TeV)')
                    axz.legend()
                    axz.text(0.98, 0.95, eff_text, transform=axz.transAxes, fontsize=8, va='top', ha='right', bbox=bbox)

                    if truth_density.sum() > 0:
                        for i, (mv, dens) in enumerate(method_density.items()):
                            ratio = np.full_like(truth_density[mask], np.nan)
                            mask2 = truth_density[mask] > 0
                            if np.any(mask2):
                                ratio[mask2] = dens[mask][mask2] / truth_density[mask][mask2]
                            axzr.step(zoom_centers, ratio, where='mid', color=colors[i % len(colors)])
                    axzr.axhline(1.0, color='k', linestyle='--', linewidth=1)
                    axzr.set_ylabel('Ratio')
                    axzr.set_xlabel('Mass [TeV]')
                    fout_zoom = os.path.join(outdir, f'mttbar_methods_vs_truth_{region}_{suffix_safe}_norm_ratio_zoom_0-{int(zoom_hi_plot*scale_mass)}.png')
                    figz.savefig(fout_zoom, dpi=150, bbox_inches='tight')
                    print('Saved zoomed normalized per-method overlay with ratio to', fout_zoom)
                    plt.close(figz)

                    # zoomed counts
                    figz2, (axz2, axz2r) = plt.subplots(2, 1, sharex=True, figsize=(8, 8), gridspec_kw={'height_ratios': [3, 1]})
                    figz2.subplots_adjust(hspace=0.05)
                    for i, (mv, cts) in enumerate(method_counts.items()):
                        cts_sub = cts[mask]
                        label = mv.replace('mttbar_', '')
                        key = mv_to_key(mv)
                        if key in eff:
                            label = f"{label} (eff={eff[key]:.3f})"
                        axz2.step(zoom_centers, cts_sub, where='mid', label=label, color=colors[i % len(colors)])
                    if truth_counts.sum() > 0:
                        truth_sub_cts = truth_counts[mask]
                        axz2.fill_between(zoom_centers, truth_sub_cts, step='mid', color='0.15', alpha=0.35, label='Truth (llbb_mass)')
                        axz2.step(zoom_centers, truth_sub_cts, where='mid', label='Truth (llbb_mass)', color='k', linewidth=1.0)
                        axz2.set_ylabel('Counts')
                        axz2.set_title(f'm(ttbar) per-method vs truth — {region} ({suffix_str}) (zoom {zoom_lo_plot:.2f}-{zoom_hi_plot:.2f} TeV)')
                    axz2.legend()
                    axz2.text(0.98, 0.95, eff_text, transform=axz2.transAxes, fontsize=8, va='top', ha='right', bbox=bbox)

                    if truth_counts.sum() > 0:
                        for i, (mv, cts) in enumerate(method_counts.items()):
                            ratio = np.full_like(truth_sub_cts, np.nan)
                            mask3 = truth_sub_cts > 0
                            if np.any(mask3):
                                ratio[mask3] = cts[mask][mask3] / truth_sub_cts[mask3]
                            axz2r.step(zoom_centers, ratio, where='mid', color=colors[i % len(colors)])
                    axz2r.axhline(1.0, color='k', linestyle='--', linewidth=1)
                    axz2r.set_ylabel('Ratio')
                    axz2r.set_xlabel('Mass [TeV]')
                    fout_zoom_counts = os.path.join(outdir, f'mttbar_methods_vs_truth_{region}_{suffix_safe}_counts_ratio_zoom_0-{int(zoom_hi_plot*scale_mass)}.png')
                    figz2.savefig(fout_zoom_counts, dpi=150, bbox_inches='tight')
                    print('Saved zoomed counts per-method overlay with ratio to', fout_zoom_counts)
                    plt.close(figz2)
        except Exception:
            pass

        # --- betaStar plots: per-method vs truth (same style as m(ttbar)) ---
        beta_method_vars = [
            'betaStar_HYPER_NUFLOWS', 'betaStar_CHI2_NUFLOWS', 'betaStar_PT_NUFLOWS', 'betaStar_DELTAR_NUFLOWS',
            'betaStar_HYPER_NUWEIGHT', 'betaStar_CHI2_NUWEIGHT', 'betaStar_PT_NUWEIGHT', 'betaStar_DELTAR_NUWEIGHT'
        ]
        legacy_beta_nuflows_map = {
            'betaStar_HYPER_NUFLOWS': 'betaStar_HYPER',
            'betaStar_CHI2_NUFLOWS': 'betaStar_CHI2',
            'betaStar_PT_NUFLOWS': 'betaStar_PT',
            'betaStar_DELTAR_NUFLOWS': 'betaStar_DELTAR',
        }
        beta_method_data = {}
        for mv in beta_method_vars:
            pats = [f'{mv}_{region}_{mult}', f'{mv}_{region}', f'{mv}_NOSYS', mv]
            legacy = legacy_beta_nuflows_map.get(mv)
            if legacy is not None:
                pats.extend([f'{legacy}_{region}_{mult}', f'{legacy}_{region}', f'{legacy}_NOSYS', legacy])
            cand = None
            for p in pats:
                if p in ks:
                    cand = p
                    break
                if (p + ';1') in ks:
                    cand = p + ';1'
                    break
            if not cand:
                continue
            obj = tree[cand]
            try:
                if hasattr(obj, 'to_numpy'):
                    counts, edges = obj.to_numpy()
                    beta_method_data[mv] = ('hist', counts, edges)
                else:
                    arr = obj.array(library='ak')
                    beta_method_data[mv] = ('ak', arr)
            except Exception:
                try:
                    arr = obj.array(library='np')
                    beta_method_data[mv] = ('array', np.asarray(arr))
                except Exception:
                    continue

        # truth for betaStar: prefer betaStar_truth_NOSYS or betaStar_truth
        beta_truth_data = None
        beta_truth_pats = [
            f'betaStar_TRUTH_{region}_{mult}', f'betaStar_TRUTH_{region}', 'betaStar_TRUTH_NOSYS', 'betaStar_TRUTH',
            f'betaStar_truth_{region}_{mult}', f'betaStar_truth_{region}', 'betaStar_truth_NOSYS', 'betaStar_truth',
            f'betaStar_{region}_{mult}', f'betaStar_{region}', 'betaStar_NOSYS', 'betaStar'
        ]
        for tpat in beta_truth_pats:
            cand = None
            if tpat in ks:
                cand = tpat
            elif (tpat + ';1') in ks:
                cand = tpat + ';1'
            if not cand:
                continue
            obj = tree[cand]
            try:
                if hasattr(obj, 'to_numpy'):
                    truth_counts, truth_edges = obj.to_numpy()
                    beta_truth_data = ('hist', truth_counts, truth_edges)
                else:
                    arr = obj.array(library='ak')
                    beta_truth_data = ('ak', arr)
                break
            except Exception:
                try:
                    tarr = np.asarray(obj.array(library='np'))
                    beta_truth_data = ('array', tarr)
                    break
                except Exception:
                    continue

        beta_method_density = {}
        beta_truth_density = np.array([])
        beta_bin_centers = np.array([])
        beta_method_counts = {}
        beta_truth_counts = np.array([])

        if not beta_method_data and beta_truth_data is None:
            print('No betaStar per-method or truth branches available for plotting for region', region, 'multiplicity', mult)
        else:
            lows = []
            highs = []
            for info in beta_method_data.values():
                if info[0] == 'hist':
                    _, edges = info[1], info[2]
                    lows.append(edges[0]); highs.append(edges[-1])
                elif info[0] == 'ak':
                    arr = info[1]
                    try:
                        npv = ak.to_numpy(arr)
                        if npv.size > 0:
                            lows.append(float(np.nanmin(npv))); highs.append(float(np.nanmax(npv)))
                    except Exception:
                        continue
                else:
                    arr = info[1]
                    try:
                        if arr.size > 0:
                            lows.append(float(np.nanmin(arr))); highs.append(float(np.nanmax(arr)))
                    except Exception:
                        continue
            if beta_truth_data is not None:
                if beta_truth_data[0] == 'hist':
                    lows.append(beta_truth_data[2][0]); highs.append(beta_truth_data[2][-1])
                elif beta_truth_data[0] == 'ak':
                    try:
                        npv = ak.to_numpy(beta_truth_data[1])
                        if npv.size > 0:
                            lows.append(float(np.nanmin(npv))); highs.append(float(np.nanmax(npv)))
                    except Exception:
                        pass
                else:
                    arr = beta_truth_data[1]
                    try:
                        if arr.size > 0:
                            lows.append(float(np.nanmin(arr))); highs.append(float(np.nanmax(arr)))
                    except Exception:
                        pass
            lo = min(lows) if lows else 0.0
            hi = max(highs) if highs else 1.0
            if lo == hi:
                lo = 0.0
                hi = max(1.0, hi)
            target_edges_beta = np.linspace(lo, hi, args.nbins + 1)

            beta_method_counts = {}
            for mv, info in beta_method_data.items():
                if info[0] == 'hist':
                    counts, edges = info[1], info[2]
                    beta_method_counts[mv] = rebin_to_target(counts, edges, target_edges_beta)
                elif info[0] == 'ak':
                    try:
                        arr = ak.to_numpy(info[1])
                        if njets_arr is not None:
                            mask = multiplicity_mask(njets_arr, mult)
                            arr = arr[mask]
                        cts, _ = np.histogram(arr, bins=target_edges_beta)
                        beta_method_counts[mv] = cts
                    except Exception:
                        beta_method_counts[mv] = np.zeros(len(target_edges_beta) - 1)
                else:
                    try:
                        arr = np.asarray(info[1])
                        cts, _ = np.histogram(arr, bins=target_edges_beta)
                        beta_method_counts[mv] = cts
                    except Exception:
                        beta_method_counts[mv] = np.zeros(len(target_edges_beta) - 1)

            if beta_truth_data is not None:
                if beta_truth_data[0] == 'hist':
                    beta_truth_counts = rebin_to_target(beta_truth_data[1], beta_truth_data[2], target_edges_beta)
                elif beta_truth_data[0] == 'ak':
                    try:
                        tarr = ak.to_numpy(beta_truth_data[1])
                        if njets_arr is not None:
                            mask = multiplicity_mask(njets_arr, mult)
                            tarr = tarr[mask]
                        beta_truth_counts, _ = np.histogram(tarr, bins=target_edges_beta)
                    except Exception:
                        beta_truth_counts = np.zeros(len(target_edges_beta) - 1)
                else:
                    try:
                        beta_truth_counts, _ = np.histogram(beta_truth_data[1], bins=target_edges_beta)
                    except Exception:
                        beta_truth_counts = np.zeros(len(target_edges_beta) - 1)
            else:
                beta_truth_counts = np.zeros(len(target_edges_beta) - 1)

            beta_bin_centers = 0.5 * (target_edges_beta[:-1] + target_edges_beta[1:])
            beta_method_density = {mv: (cts / np.sum(cts) if np.sum(cts) > 0 else cts) for mv, cts in beta_method_counts.items()}
            beta_truth_density = beta_truth_counts / np.sum(beta_truth_counts) if np.sum(beta_truth_counts) > 0 else beta_truth_counts

            fig, (ax, axr) = plt.subplots(2, 1, sharex=True, figsize=(8, 8), gridspec_kw={'height_ratios': [3, 1]})
            fig.subplots_adjust(hspace=0.05)
            for i, (mv, dens) in enumerate(beta_method_density.items()):
                key = mv_to_key(mv)
                label = mv.replace('betaStar_', '')
                if key in eff:
                    label = f"{label} (eff={eff[key]:.3f})"
                ax.step(beta_bin_centers, dens, where='mid', label=label, color=colors[i % len(colors)])
            if beta_truth_counts.sum() > 0:
                ax.fill_between(beta_bin_centers, beta_truth_density, step='mid', color='0.15', alpha=0.35, label='Truth (betaStar)')
                ax.step(beta_bin_centers, beta_truth_density, where='mid', label='Truth (betaStar)', color='k', linewidth=1.0)
            ax.set_ylabel('Normalized density')
            ax.set_title(f'betaStar per-method vs truth — {region} ({suffix_str}) (normalized)')
            ax.legend()
            ax.text(0.98, 0.95, eff_text, transform=ax.transAxes, fontsize=8, va='top', ha='right', bbox=bbox)

            if beta_truth_density.sum() > 0:
                for i, (mv, dens) in enumerate(beta_method_density.items()):
                    ratio = np.full_like(beta_truth_density, np.nan)
                    mask = beta_truth_density > 0
                    if np.any(mask):
                        ratio[mask] = dens[mask] / beta_truth_density[mask]
                    axr.step(beta_bin_centers, ratio, where='mid', color=colors[i % len(colors)])
                axr.axhspan(0.8, 1.2, color='gray', alpha=0.12)
                axr.set_ylim(0.5, 1.5)
            axr.axhline(1.0, color='k', linestyle='--', linewidth=1)
            axr.set_ylabel('Ratio')
            axr.set_xlabel('betaStar')

            fout_norm_ratio = os.path.join(outdir, f'betastar_methods_vs_truth_{region}_{suffix_safe}_norm_ratio.png')
            fig.savefig(fout_norm_ratio, dpi=150, bbox_inches='tight')
            print('Saved normalized per-method betaStar overlay with ratio to', fout_norm_ratio)
            plt.close(fig)

            fig2, (ax2, ax2r) = plt.subplots(2, 1, sharex=True, figsize=(8, 8), gridspec_kw={'height_ratios': [3, 1]})
            fig2.subplots_adjust(hspace=0.05)
            for i, (mv, cts) in enumerate(beta_method_counts.items()):
                key = mv_to_key(mv)
                label = mv.replace('betaStar_', '')
                if key in eff:
                    label = f"{label} (eff={eff[key]:.3f})"
                ax2.step(beta_bin_centers, cts, where='mid', label=label, color=colors[i % len(colors)])
            if beta_truth_counts.sum() > 0:
                ax2.fill_between(beta_bin_centers, beta_truth_counts, step='mid', color='0.15', alpha=0.35, label='Truth (betaStar)')
                ax2.step(beta_bin_centers, beta_truth_counts, where='mid', label='Truth (betaStar)', color='k', linewidth=1.0)
            ax2.set_ylabel('Counts')
            ax2.set_title(f'betaStar per-method vs truth — {region} ({suffix_str}) (counts)')
            ax2.legend()
            ax2.text(0.98, 0.95, eff_text, transform=ax2.transAxes, fontsize=8, va='top', ha='right', bbox=bbox)

            if beta_truth_counts.sum() > 0:
                for i, (mv, cts) in enumerate(beta_method_counts.items()):
                    ratio = np.full_like(beta_truth_counts, np.nan)
                    mask = beta_truth_counts > 0
                    if np.any(mask):
                        ratio[mask] = cts[mask] / beta_truth_counts[mask]
                    ax2r.step(beta_bin_centers, ratio, where='mid', color=colors[i % len(colors)])
                ax2r.axhspan(0.8, 1.2, color='gray', alpha=0.12)
                ax2r.set_ylim(0.5, 1.5)
            ax2r.axhline(1.0, color='k', linestyle='--', linewidth=1)
            ax2r.set_ylabel('Ratio')
            ax2r.set_xlabel('betaStar')

            fout_counts_ratio = os.path.join(outdir, f'betastar_methods_vs_truth_{region}_{suffix_safe}_counts_ratio.png')
            fig2.savefig(fout_counts_ratio, dpi=150, bbox_inches='tight')
            print('Saved counts per-method betaStar overlay with ratio to', fout_counts_ratio)
            plt.close(fig2)

        # return collected efficiencies/totals and per-multiplicity histograms/densities
        return (eff, err, totals,
            method_density, truth_density, bin_centers, method_counts, truth_counts,
            beta_method_density, beta_truth_density, beta_bin_centers, beta_method_counts, beta_truth_counts,
            suffix_safe)

    # Main: iterate requested multiplicities and collect efficiencies
    jet_mults = [2, 3, 4, 5, '6+']
    any_hist_found = False
    results_by_mult = {}

    # Also keep an overall (inclusive) view for superplots.
    overall_res = process_multiplicity('')
    if overall_res:
        results_by_mult['overall'] = overall_res

    for m in jet_mults:
        suf = find_suffix_for_mult(m)
        if suf:
            any_hist_found = True
            res = process_multiplicity(suf)
        else:
            print('No per-multiplicity histogram found for', m, '- attempting ntuple/event-level processing')
            res = process_multiplicity(m)
        if res:
            # store results under canonical multiplicity label
            results_by_mult[str(m)] = res
            any_hist_found = True

    # If we collected efficiencies across multiplicities, make an efficiency vs multiplicity plot
    if results_by_mult:
        try:
            import matplotlib
            matplotlib.use('Agg')
            import matplotlib.pyplot as plt
            mult_labels = [str(x) for x in jet_mults]
            x = np.arange(len(mult_labels))
            # collect method names
            method_set = set()
            for mres in results_by_mult.values():
                effd = mres[0]
                method_set.update(list(effd.keys()))
            methods = sorted(method_set)
            fig, ax = plt.subplots(figsize=(7,4))
            for i, meth in enumerate(methods):
                ys = []
                yerr = []
                for ml in mult_labels:
                    if ml in results_by_mult and meth in results_by_mult[ml][0]:
                        ys.append(results_by_mult[ml][0][meth])
                        yerr.append(results_by_mult[ml][1][meth])
                    else:
                        ys.append(np.nan)
                        yerr.append(0.0)
                ax.errorbar(x, ys, yerr=yerr, marker='o', label=meth)
            ax.set_xticks(x)
            ax.set_xticklabels(mult_labels)
            ax.set_ylim(0,1)
            ax.set_xlabel('Jet multiplicity')
            ax.set_ylabel('Reconstruction efficiency')
            ax.set_title(f'Efficiency vs jet multiplicity — {region}')
            ax.legend()
            fout_ev = os.path.join(outdir, f'efficiency_vs_multiplicity_{region}.png')
            fig.savefig(fout_ev, dpi=150, bbox_inches='tight')
            print('Saved efficiency-vs-multiplicity plot to', fout_ev)
            # --- Superplots: include overall + multiplicities ---
            try:
                import matplotlib.pyplot as plt

                mult_order = ['overall', '2', '3', '4', '5', '6+']
                colors_sp = ['#1f77b4','#ff7f0e','#2ca02c','#d62728']
                style_map_sp = {'baseline': '-', 'chisquared': '--', 'deltaR': '-.', 'hyper': ':'}

                def mv_to_key_local(mv):
                    l = mv.lower()
                    if 'chi' in l: return 'chisquared'
                    if 'delt' in l or 'delta' in l: return 'deltaR'
                    if 'pt' in l: return 'baseline'
                    if 'hyper' in l: return 'hyper'
                    return None

                # ---------- m(ttbar) SUPERPLOT (normalized) ----------
                fig, axes = plt.subplots(2, 3, figsize=(16, 10), sharex=True, sharey=True)
                axes = axes.flatten()

                for idx, ml in enumerate(mult_order):
                    ax = axes[idx]

                    if ml not in results_by_mult:
                        ax.set_visible(False)
                        continue

                    res = results_by_mult[ml]
                    eff, err = res[0], res[1]
                    method_density = res[3]
                    truth_density = res[4]
                    bin_centers = res[5]

                    # efficiency text
                    eff_text = "\n".join([f"{k}: {eff[k]:.3f}" for k in eff])

                    for i, (mv, dens) in enumerate(method_density.items()):
                        key = mv_to_key_local(mv)
                        ls = style_map_sp.get(key, '-')
                        label = mv.replace('mttbar_', '')
                        ax.step(bin_centers, dens, where='mid',
                                color=colors_sp[i % len(colors_sp)],
                                linestyle=ls,
                                label=label)

                    if truth_density.sum() > 0:
                        ax.fill_between(bin_centers, truth_density, step='mid',
                                        color='0.15', alpha=0.3)
                        ax.step(bin_centers, truth_density, where='mid', color='k')

                    ax.set_title(f'{ml} jets')
                    ax.text(0.97, 0.97, eff_text, transform=ax.transAxes,
                            ha='right', va='top', fontsize=8,
                            bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))

                    ax.grid(True, linestyle=':', alpha=0.5)

                fig.suptitle(f'm(ttbar) comparison — {region} (all multiplicities)')
                fig.legend(loc='upper right')

                fout = os.path.join(outdir, f'mttbar_superplot_all_{region}.png')
                fig.savefig(fout, dpi=150, bbox_inches='tight')
                plt.close(fig)

                print('Saved mttbar superplot:', fout)

                # ---------- betaStar SUPERPLOT ----------
                fig, axes = plt.subplots(2, 3, figsize=(16, 10), sharex=True, sharey=True)
                axes = axes.flatten()

                for idx, ml in enumerate(mult_order):
                    ax = axes[idx]

                    if ml not in results_by_mult:
                        ax.set_visible(False)
                        continue

                    res = results_by_mult[ml]
                    eff = res[0]

                    beta_method_density = res[8]
                    beta_truth_density = res[9]
                    beta_bin_centers = res[10]

                    for i, (mv, dens) in enumerate(beta_method_density.items()):
                        key = mv_to_key_local(mv)
                        ls = style_map_sp.get(key, '-')
                        label = mv.replace('betaStar_', '')
                        ax.step(beta_bin_centers, dens, where='mid',
                                color=colors_sp[i % len(colors_sp)],
                                linestyle=ls,
                                label=label)

                    if beta_truth_density.size > 0 and np.nansum(beta_truth_density) > 0:
                        ax.fill_between(beta_bin_centers, beta_truth_density, step='mid',
                                        color='0.15', alpha=0.3)
                        ax.step(beta_bin_centers, beta_truth_density, where='mid', color='k')

                    ax.set_title(f'{ml} jets')
                    ax.text(0.97, 0.97,
                            "\n".join([f"{k}: {eff[k]:.3f}" for k in eff]),
                            transform=ax.transAxes,
                            ha='right', va='top', fontsize=8,
                            bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))

                    ax.grid(True, linestyle=':', alpha=0.5)

                fig.suptitle(f'betaStar comparison — {region} (all multiplicities)')
                fout = os.path.join(outdir, f'betastar_superplot_all_{region}.png')
                fig.savefig(fout, dpi=150, bbox_inches='tight')
                plt.close(fig)

                print('Saved betaStar superplot:', fout)

            except Exception as e:
                print('Superplot generation failed:', e)
        except Exception:
            pass

    # If nothing was found at all, try a full-region fallback (no multiplicity)
    if not any_hist_found:
        # check for at least one of the full-region keys
        def has_full_region_keys():
            candidates = [f'mttbar_HYPER_{region}', f'llbb_mass_{region}', f'baseline_pairing_correct_{region}', 'invariant_mass_llbb_NOSYS', 'mttbar_HYPER_NOSYS']
            for c in candidates:
                if c in ks or (c + ';1') in ks:
                    return True
            return False

        if has_full_region_keys():
            print('No multiplicity-specific histograms found; producing full-region plots')
            process_multiplicity('')
        else:
            print('No full-region keys either; nothing to plot.')

if __name__ == '__main__':
    main()
