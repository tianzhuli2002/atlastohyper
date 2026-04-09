#!/usr/bin/env python3
import argparse
import os
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
import awkward as ak

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
                if k.lower().endswith(f'_{region}_{cand}'.lower()):
                    return cand
        return None

    # Colors & Linestyles
    # NUFLOWS (saturated) vs NUWEIGHT (lighter shade)
    colors_nf = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728']
    colors_nw = ['#6baed6', '#fdbb84', '#74c476', '#e88a8a']
    styles = {'baseline': '-', 'chisquared': '--', 'deltaR': '-.', 'hyper': ':'}

    def mv_to_key(mv):
        l = mv.lower()
        if 'chi' in l: return 'chisquared'
        if 'delt' in l or 'delta' in l: return 'deltaR'
        if 'pt' in l: return 'baseline'
        if 'hyper' in l: return 'hyper'
        return None

    def get_color_style(mv):
        key = mv_to_key(mv)
        if key in styles:
            idx = list(styles.keys()).index(key)
            col = colors_nw[idx] if 'nuweight' in mv.lower() else colors_nf[idx]
            return col, styles[key]
        return '#000000', '-'

    def process_multiplicity(mult):
        import math

        method_patterns = {
            'baseline': [f'PT_correct_{region}_{mult}', f'PT_correct_{region}', 'PT_correct_NOSYS', 'PT_correct'],
            'chisquared': [f'CHI2_correct_{region}_{mult}', f'CHI2_correct_{region}', 'CHI2_correct_NOSYS', 'CHI2_correct'],
            'deltaR': [f'DELTAR_correct_{region}_{mult}', f'DELTAR_correct_{region}', 'DELTAR_correct_NOSYS', 'DELTAR_correct'],
            'hyper': [f'HYPER_correct_{region}_{mult}', f'HYPER_correct_{region}', 'HYPER_correct_NOSYS', 'HYPER_correct']
        }

        found = {}
        for m, pats in method_patterns.items():
            candidate = None
            for p in pats:
                if p in ks: candidate = p; break
                if (p + ';1') in ks: candidate = p + ';1'; break
            if candidate is None: continue
            obj = tree[candidate]
            try:
                if hasattr(obj, 'to_numpy'):
                    counts, edges = obj.to_numpy()
                    found[m] = ('hist', counts, edges)
                else:
                    found[m] = ('ak', obj.array(library='ak'))
            except Exception:
                continue

        if not found:
            print(f'No method correctness branches found for multiplicity {mult}')
            return None

        def get_njets_array():
            jet_cands = ['jet_pt_NOSYS', 'jet_pt', 'Jet_pt', 'Jet_pt_NOSYS', 'jet_TLV_NOSYS', 'jet_TLV']
            for cand in jet_cands:
                k = cand if cand in ks else (cand + ';1' if (cand + ';1') in ks else None)
                if not k: continue
                try:
                    arr_ak = tree[k].array(library='ak')
                    nums = ak.num(arr_ak)
                    return np.asarray(ak.to_numpy(nums), dtype=int)
                except Exception: continue
            nj_cands = ['n_jets_NOSYS', 'n_jets', 'nJets', 'nJets_good_NOSYS', 'nJets_good']
            for cand in nj_cands:
                k = cand if cand in ks else (cand + ';1' if (cand + ';1') in ks else None)
                if not k: continue
                try: return np.asarray(tree[k].array(library='np'), dtype=int)
                except Exception: continue
            return None

        njets_arr = get_njets_array()

        def multiplicity_mask(njets, mult_label):
            if njets is None: return None
            arr = np.asarray(njets)
            s = str(mult_label)
            if s.endswith('+'):
                try: return arr >= int(s[:-1])
                except: return np.ones(arr.shape[0], dtype=bool)
            try:
                n = int(s.lstrip('_').split('j')[0])
                return arr == n
            except: return np.ones(arr.shape[0], dtype=bool)

        eff, err, totals = {}, {}, {}
        for m, info in found.items():
            if info[0] == 'hist':
                counts, edges = info[1], info[2]
                centers = 0.5 * (edges[:-1] + edges[1:])
                idx_true = int(np.argmin(np.abs(centers - 1.0)))
                n_true = float(counts[idx_true])
                total = float(np.sum(counts))
            else:
                try:
                    vals_all = ak.to_numpy(info[1])
                    if njets_arr is not None:
                        mask = multiplicity_mask(njets_arr, mult)
                        vals = np.asarray(vals_all)[mask]
                    else:
                        vals = np.asarray(vals_all)
                    valid_mask = (vals == 0) | (vals == 1)
                    vals = vals[valid_mask]
                    total = float(vals.size)
                    n_true = float(np.count_nonzero(vals == 1)) if total > 0 else 0.0
                except Exception:
                    n_true = total = 0.0

            e = n_true / total if total > 0 else 0.0
            ee = math.sqrt(e * (1.0 - e) / total) if total > 0 else 0.0
            eff[m] = e
            err[m] = ee
            totals[m] = int(total)

        print(f'Multiplicity {mult}: {", ".join(f"{m}: {eff[m]:.3f}" for m in eff)}')

        # Generic function to collect and process method/truth data
        def collect_data(prefix, method_vars, legacy_map):
            method_data = {}
            for mv in method_vars:
                pats = [f'{mv}_{region}_{mult}', f'{mv}_{region}', f'{mv}_NOSYS', mv]
                lg = legacy_map.get(mv)
                if lg: pats.extend([f'{lg}_{region}_{mult}', f'{lg}_{region}', f'{lg}_NOSYS', lg])
                cand = next((p for p in pats if p in ks or (p+';1') in ks), None)
                if not cand: continue
                cand = cand if cand in ks else cand + ';1'
                obj = tree[cand]
                try:
                    if hasattr(obj, 'to_numpy'): method_data[mv] = ('hist', *obj.to_numpy())
                    else: method_data[mv] = ('ak', obj.array(library='ak'))
                except Exception: continue

            # Truth data
            truth_data = None
            for tpat in [f'{prefix}_TRUTH_{region}_{mult}', f'{prefix}_TRUTH_{region}', f'{prefix}_TRUTH_NOSYS', f'{prefix}_TRUTH',
                         f'{prefix}_truth_{region}_{mult}', f'{prefix}_truth_{region}', f'{prefix}_truth_NOSYS', f'{prefix}_truth']:
                cand = tpat if tpat in ks else (tpat+';1' if (tpat+';1') in ks else None)
                if not cand: continue
                obj = tree[cand]
                try:
                    if hasattr(obj, 'to_numpy'): truth_data = ('hist', *obj.to_numpy())
                    else: truth_data = ('ak', obj.array(library='ak'))
                    break
                except Exception: continue

            if not method_data: return None

            # Unified binning
            lows, highs = [], []
            for info in method_data.values():
                if info[0] == 'hist': lows.append(info[2][0]); highs.append(info[2][-1])
                elif info[0] == 'ak':
                    try:
                        v = ak.to_numpy(info[1]); lows.append(np.nanmin(v)); highs.append(np.nanmax(v))
                    except: continue
            if truth_data:
                if truth_data[0] == 'hist': lows.append(truth_data[2][0]); highs.append(truth_data[2][-1])
                elif truth_data[0] == 'ak':
                    try: v = ak.to_numpy(truth_data[1]); lows.append(np.nanmin(v)); highs.append(np.nanmax(v))
                    except: pass
            
            lo, hi = min(lows), max(highs)
            # For mttbar, cap at args.xmax (GeV) * 1000 to match MeV inputs
            if prefix == 'mttbar': hi = min(hi, args.xmax * 1000.0)
            if lo == hi: lo, hi = 0.0, 1.0
            target_edges = np.linspace(lo, hi, args.nbins + 1)

            # Build counts
            method_counts = {}
            for mv, info in method_data.items():
                if info[0] == 'hist':
                    method_counts[mv] = rebin_to_target(info[1], info[2], target_edges)
                elif info[0] == 'ak':
                    arr = ak.to_numpy(info[1])
                    if njets_arr is not None: arr = arr[multiplicity_mask(njets_arr, mult)]
                    method_counts[mv] = np.histogram(arr, bins=target_edges)[0]

            if truth_data:
                if truth_data[0] == 'hist':
                    truth_counts = rebin_to_target(truth_data[1], truth_data[2], target_edges)
                elif truth_data[0] == 'ak':
                    tarr = ak.to_numpy(truth_data[1])
                    if njets_arr is not None: tarr = tarr[multiplicity_mask(njets_arr, mult)]
                    truth_counts = np.histogram(tarr, bins=target_edges)[0]
            else:
                truth_counts = np.zeros(len(target_edges)-1)

            # Scale mass to GeV if mttbar
            if prefix in ('mttbar', 'llbb'):
                scale = 1000.0
                bin_centers = 0.5 * (target_edges[:-1] + target_edges[1:]) / scale
            else:
                scale = 1.0
                bin_centers = 0.5 * (target_edges[:-1] + target_edges[1:])
            
            method_density = {mv: c/np.sum(c) if np.sum(c)>0 else c for mv, c in method_counts.items()}
            truth_density = truth_counts/np.sum(truth_counts) if np.sum(truth_counts)>0 else truth_counts

            return {
                'density': method_density, 'truth': truth_density, 
                'centers': bin_centers, 'counts': method_counts, 
                'truth_counts': truth_counts, 'edges': target_edges
            }

        def collect_triplet_data(method_name_1, method_name_2, truth_name):
            def read_branch(branch_name):
                pats = [f'{branch_name}_{region}_{mult}', f'{branch_name}_{region}', f'{branch_name}_NOSYS', branch_name]
                cand = next((p for p in pats if p in ks or (p + ';1') in ks), None)
                if not cand:
                    return None
                cand = cand if cand in ks else cand + ';1'
                obj = tree[cand]
                try:
                    if hasattr(obj, 'to_numpy'):
                        return ('hist', *obj.to_numpy())
                    return ('ak', obj.array(library='ak'))
                except Exception:
                    return None

            d1 = read_branch(method_name_1)
            d2 = read_branch(method_name_2)
            dt = read_branch(truth_name)
            if d1 is None and d2 is None and dt is None:
                return None

            all_data = []
            for d in [d1, d2, dt]:
                if d is not None:
                    all_data.append(d)

            lows, highs = [], []
            for info in all_data:
                if info[0] == 'hist':
                    lows.append(info[2][0])
                    highs.append(info[2][-1])
                else:
                    try:
                        arr = ak.to_numpy(info[1])
                        if njets_arr is not None:
                            arr = arr[multiplicity_mask(njets_arr, mult)]
                        arr = arr[np.isfinite(arr)]
                        if arr.size > 0:
                            lows.append(np.min(arr))
                            highs.append(np.max(arr))
                    except Exception:
                        continue
            if not lows or not highs:
                return None

            lo, hi = min(lows), max(highs)
            if lo == hi:
                lo, hi = 0.0, 1.0
            target_edges = np.linspace(lo, hi, args.nbins + 1)

            def to_counts(info):
                if info is None:
                    return np.zeros(len(target_edges) - 1)
                if info[0] == 'hist':
                    return rebin_to_target(info[1], info[2], target_edges)
                arr = ak.to_numpy(info[1])
                if njets_arr is not None:
                    arr = arr[multiplicity_mask(njets_arr, mult)]
                arr = arr[np.isfinite(arr)]
                return np.histogram(arr, bins=target_edges)[0]

            c1 = to_counts(d1)
            c2 = to_counts(d2)
            ct = to_counts(dt)

            return {
                'density': {
                    method_name_1.upper(): c1 / np.sum(c1) if np.sum(c1) > 0 else c1,
                    method_name_2.upper(): c2 / np.sum(c2) if np.sum(c2) > 0 else c2,
                },
                'truth': ct / np.sum(ct) if np.sum(ct) > 0 else ct,
                'centers': 0.5 * (target_edges[:-1] + target_edges[1:]),
                'counts': {
                    method_name_1.upper(): c1,
                    method_name_2.upper(): c2,
                },
                'truth_counts': ct,
                'edges': target_edges,
            }

        # mttbar
        mtt_vars = ['mttbar_HYPER_NUFLOWS', 'mttbar_CHI2_NUFLOWS', 'mttbar_PT_NUFLOWS', 'mttbar_DELTAR_NUFLOWS',
                    'mttbar_HYPER_NUWEIGHT', 'mttbar_CHI2_NUWEIGHT', 'mttbar_PT_NUWEIGHT', 'mttbar_DELTAR_NUWEIGHT']
        mtt_legacy = {'mttbar_HYPER_NUFLOWS': 'mttbar_HYPER', 'mttbar_CHI2_NUFLOWS': 'mttbar_CHI2',
                      'mttbar_PT_NUFLOWS': 'mttbar_PT', 'mttbar_DELTAR_NUFLOWS': 'mttbar_DELTAR'}
        mtt_data = collect_data('mttbar', mtt_vars, mtt_legacy)

        # betaStar
        beta_vars = ['betaStar_HYPER_NUFLOWS', 'betaStar_CHI2_NUFLOWS', 'betaStar_PT_NUFLOWS', 'betaStar_DELTAR_NUFLOWS',
                     'betaStar_HYPER_NUWEIGHT', 'betaStar_CHI2_NUWEIGHT', 'betaStar_PT_NUWEIGHT', 'betaStar_DELTAR_NUWEIGHT']
        beta_legacy = {'betaStar_HYPER_NUFLOWS': 'betaStar_HYPER', 'betaStar_CHI2_NUFLOWS': 'betaStar_CHI2',
                       'betaStar_PT_NUFLOWS': 'betaStar_PT', 'betaStar_DELTAR_NUFLOWS': 'betaStar_DELTAR'}
        beta_data = collect_data('betaStar', beta_vars, beta_legacy)

        # llbb
        llbb_vars = ['llbb_HYPER_NUFLOWS', 'llbb_CHI2_NUFLOWS', 'llbb_PT_NUFLOWS', 'llbb_DELTAR_NUFLOWS',
                     'llbb_HYPER_NUWEIGHT', 'llbb_CHI2_NUWEIGHT', 'llbb_PT_NUWEIGHT', 'llbb_DELTAR_NUWEIGHT']
        llbb_legacy = {}
        llbb_data = collect_data('llbb', llbb_vars, llbb_legacy)

        # neutrino kinematics: NUFLOWS and NUWEIGHT compared to truth
        nu1_data = collect_triplet_data('nuflows_nu1_pt', 'nu_pt', 'truth_nu1_pt')
        nu2_data = collect_triplet_data('nuflows_nu2_pt', 'nubar_pt', 'truth_nu2_pt')

        return {
            'eff': eff, 'err': err, 'totals': totals,
            'mttbar': mtt_data, 'betaStar': beta_data,
            'llbb': llbb_data,
            'nu1_pt': nu1_data,
            'nu2_pt': nu2_data,
        }

    # Collect data across multiplicities
    mult_order = ['overall', '2', '3', '4', '5', '6+']
    results_by_mult = {}

    # Audit mttbar availability and fill quality before plotting.
    mttbar_audit_branches = [
        'mttbar_HYPER_NUFLOWS_NOSYS',
        'mttbar_HYPER_NUWEIGHT_NOSYS',
        'mttbar_CHI2_NUFLOWS_NOSYS',
        'mttbar_CHI2_NUWEIGHT_NOSYS',
        'mttbar_PT_NUFLOWS_NOSYS',
        'mttbar_PT_NUWEIGHT_NOSYS',
        'mttbar_DELTAR_NUFLOWS_NOSYS',
        'mttbar_DELTAR_NUWEIGHT_NOSYS',
        'mttbar_TRUTH_NOSYS',
    ]
    print('mttbar branch audit:')
    for b in mttbar_audit_branches:
        bk = b if b in ks else (b + ';1' if (b + ';1') in ks else None)
        if bk is None:
            print(f'  MISSING: {b}')
            continue
        try:
            vals = np.asarray(tree[bk].array(library='np'), dtype=float)
            finite = np.isfinite(vals)
            valid = finite & (vals > 0)
            n_tot = vals.size
            n_valid = int(np.count_nonzero(valid))
            frac = (n_valid / n_tot) if n_tot > 0 else 0.0
            print(f'  {b}: present, valid(>0)={n_valid}/{n_tot} ({frac:.3f})')
        except Exception as exc:
            print(f'  {b}: present but unreadable ({exc})')

    # Process overall
    res = process_multiplicity('')
    if res: results_by_mult['overall'] = res

    # Process specific multiplicities
    for m in [2, 3, 4, 5, '6+']:
        suf = find_suffix_for_mult(m)
        if not suf: suf = str(m)
        r = process_multiplicity(suf)
        if r: results_by_mult[str(m)] = r

    if not results_by_mult:
        print('No data found for any multiplicity.')
        return

    # 1. COMBINED EFFICIENCY BAR PLOT
    methods = ['baseline', 'chisquared', 'deltaR', 'hyper']
    fig, ax = plt.subplots(figsize=(8, 5))
    x = np.arange(len(mult_order))
    width = 0.15
    for i, meth in enumerate(methods):
        heights = []
        yerrs = []
        for ml in mult_order:
            if ml in results_by_mult:
                d = results_by_mult[ml]
                h = d['eff'].get(meth, np.nan)
                e = d['err'].get(meth, 0)
                heights.append(h)
                yerrs.append(e)
            else:
                heights.append(np.nan)
                yerrs.append(0)
        ax.bar(x + i*width, heights, width, label=meth.capitalize(), color=colors_nf[i], edgecolor='black', yerr=yerrs, capsize=3)
    ax.set_xticks(x + 1.5*width)
    ax.set_xticklabels([m.replace('6+','6+') for m in mult_order])
    ax.set_ylabel('Reconstruction Efficiency')
    ax.set_xlabel('Jet Multiplicity')
    ax.set_title('Reconstruction Efficiency by Method and Multiplicity')
    ax.set_ylim(0, 1.1)
    ax.legend(loc='lower right')
    ax.grid(axis='y', linestyle=':', alpha=0.5)
    fig.tight_layout()
    fig.savefig(os.path.join(outdir, 'efficiency_combined.png'), dpi=150)
    print('Saved combined efficiency plot')
    plt.close(fig)

    # Helper for superplot generation
    def make_superplot(data_key, save_name, x_label):
        fig = plt.figure(figsize=(16, 12))
        gs = GridSpec(4, 3, figure=fig, hspace=0.05, wspace=0.08, height_ratios=[3, 1, 3, 1])

        for idx, ml in enumerate(mult_order):
            row_main = 2 * (idx // 3)
            row_ratio = row_main + 1
            col = idx % 3

            if ml not in results_by_mult or results_by_mult[ml].get(data_key) is None:
                fig.add_subplot(gs[row_main, col]).set_visible(False)
                fig.add_subplot(gs[row_ratio, col]).set_visible(False)
                continue

            d = results_by_mult[ml][data_key]
            md, td, bc = d['density'], d['truth'], d['centers']
            eff = results_by_mult[ml]['eff']

            ax_main = fig.add_subplot(gs[row_main, col])
            ax_ratio = fig.add_subplot(gs[row_ratio, col], sharex=ax_main)

            for mv, dens in md.items():
                if 'NUFLOWS_NU' in mv:
                    col_s, ls = '#1f77b4', '-'
                elif 'NU_PT' in mv or 'NUBAR_PT' in mv:
                    col_s, ls = '#d62728', '--'
                else:
                    col_s, ls = get_color_style(mv)
                ax_main.step(bc, dens, where='mid', color=col_s, linestyle=ls, label=mv.replace(f'{data_key}_',''))

            if np.sum(td) > 0:
                ax_main.fill_between(bc, td, step='mid', color='0.15', alpha=0.3)
                ax_main.step(bc, td, where='mid', color='k', linewidth=1.0)

            ax_main.set_title(f'{ml} jets', pad=8)
            ax_main.text(0.02, 0.95, "\n".join(f"{k}: {eff[k]:.3f}" for k in eff),
                         transform=ax_main.transAxes, fontsize=7, va='top', ha='left',
                         bbox=dict(boxstyle='round', facecolor='white', alpha=0.85, edgecolor='gray'))
            ax_main.grid(True, linestyle=':', alpha=0.5)
            ax_main.tick_params(labelbottom=False)

            # Legend only on first plot
            if idx == 0:
                ax_main.legend(loc='upper right', fontsize=7, framealpha=0.9)

            # Ratio panel
            mask_t = td > 0
            if np.any(mask_t):
                for mv, dens in md.items():
                    if 'NUFLOWS_NU' in mv:
                        col_s, ls = '#1f77b4', '-'
                    elif 'NU_PT' in mv or 'NUBAR_PT' in mv:
                        col_s, ls = '#d62728', '--'
                    else:
                        col_s, ls = get_color_style(mv)
                    r = np.full_like(td, np.nan)
                    r[mask_t] = dens[mask_t] / td[mask_t]
                    ax_ratio.step(bc, r, where='mid', color=col_s, linestyle=ls)
            ax_ratio.axhline(1.0, color='k', linestyle='--', linewidth=1)
            ax_ratio.axhspan(0.8, 1.2, color='gray', alpha=0.1)
            ax_ratio.set_ylim(0.5, 1.5)
            ax_ratio.grid(True, linestyle=':', alpha=0.5)
            ax_ratio.tick_params(labelsize=7)

            if row_ratio == 3:
                ax_ratio.set_xlabel(x_label, fontsize=9)
            else:
                ax_ratio.tick_params(labelbottom=False)

        # Adjust title spacing
        fig.subplots_adjust(top=0.93)
        fig.suptitle(f'{data_key} Comparison — All Multiplicities', fontsize=14, y=0.99)
        fig.savefig(os.path.join(outdir, save_name), dpi=150, bbox_inches='tight')
        print(f'Saved {save_name}')
        plt.close(fig)

    # 2. SUPERPLOTS
    make_superplot('mttbar', 'mttbar_superplot.png', 'Mass [GeV]')
    make_superplot('betaStar', 'betastar_superplot.png', 'betaStar')
    make_superplot('llbb', 'llbb_superplot.png', 'Mass [GeV]')
    make_superplot('nu1_pt', 'nu1_pt_superplot.png', 'pT [MeV]')
    make_superplot('nu2_pt', 'nu2_pt_superplot.png', 'pT [MeV]')

    # 3. ZOOM PLOTS (mttbar only, fixed scaling & visibility)
    zoom_lo, zoom_hi = 0.0, float(args.xmax)
    print(f'Generating zoom plots for range {zoom_lo}-{zoom_hi} GeV')
    
    # Iterate over multiplicities that have data
    zoom_mults = [m for m in ['2', '3', '4', '5', '6+', 'overall'] if m in results_by_mult and results_by_mult[m].get('mttbar')]
    
    for ml in zoom_mults:
        d = results_by_mult[ml]['mttbar']
        bc = d['centers']
        md, td = d['density'], d['truth']
        mask = (bc >= zoom_lo) & (bc <= zoom_hi)
        
        if not np.any(mask): 
            print(f'  Skipping {ml}: No data in zoom range')
            continue

        fig, (ax, axr) = plt.subplots(2, 1, figsize=(9, 8), sharex=True,
                                       gridspec_kw={'height_ratios': [3, 1]},
                                       subplot_kw={'facecolor': 'white'})
        fig.subplots_adjust(hspace=0.05, top=0.90)

        for mv, dens in md.items():
            col_s, ls = get_color_style(mv)
            sub = dens[mask]
            ax.step(bc[mask], sub, where='mid', color=col_s, linestyle=ls, label=mv.replace('mttbar_',''))

        if np.sum(td) > 0:
            tsub = td[mask]
            ax.fill_between(bc[mask], tsub, step='mid', color='0.15', alpha=0.35)
            ax.step(bc[mask], tsub, where='mid', color='k', linewidth=1.0)

        # Ensure y-axis limits are set so plots are visible
        max_dens = max((np.max(dens[mask]) for dens in md.values() if np.any(mask)), default=0)
        max_truth = np.max(td[mask]) if np.sum(td[mask]) > 0 else 0
        y_max = max(max_dens, max_truth)
        if y_max == 0: y_max = 0.01
        ax.set_ylim(0, y_max * 1.2)
        
        ax.set_ylabel('Normalized Density', fontsize=10)
        ax.set_title(f'm(ttbar) — {ml} jets (Zoom: {zoom_lo:.0f}-{zoom_hi:.0f} GeV)', fontsize=12)
        ax.legend(loc='upper right', fontsize=9)
        ax.grid(True, linestyle=':', alpha=0.5)
        ax.tick_params(axis='both', which='major', labelsize=9)

        # Ratio
        mask_t = td[mask] > 0
        if np.any(mask_t):
            for mv, dens in md.items():
                col_s, ls = get_color_style(mv)
                r = np.full_like(td[mask], np.nan)
                r[mask_t] = dens[mask][mask_t] / td[mask][mask_t]
                axr.step(bc[mask], r, where='mid', color=col_s, linestyle=ls)
        axr.axhline(1.0, color='k', linestyle='--', linewidth=1)
        axr.axhspan(0.8, 1.2, color='gray', alpha=0.1)
        axr.set_ylim(0.5, 1.5)
        axr.set_ylabel('Ratio', fontsize=10)
        axr.set_xlabel('Mass [GeV]', fontsize=10)
        axr.grid(True, linestyle=':', alpha=0.5)
        axr.tick_params(axis='both', which='major', labelsize=9)

        fname = f'mttbar_zoom_{ml.replace("+","plus")}jets.png'
        fig.savefig(os.path.join(outdir, fname), dpi=150, bbox_inches='tight')
        print(f'Saved zoom plot: {fname}')
        plt.close(fig)

if __name__ == '__main__':
    main()