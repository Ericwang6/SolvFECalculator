import os, shutil
import json
import subprocess
from pathlib import Path
from typing import List


def split_top(top: os.PathLike, atp: os.PathLike, itp: os.PathLike, posre: os.PathLike):
    atpf = open(atp, 'w')
    itpf = open(itp, 'w')
    with open(top, 'r') as f:
        read_atp = False
        read_itp = False
        read_atoms = False
        atoms = []
        for line in f:
            if line.startswith("[ atomtypes ]"):
                read_atp = True
                read_itp = False
            elif line.startswith("[ moleculetype ]"):
                read_atp = False
                read_itp = False
                itpf.write("[ moleculetype ]\n MOL   3\n\n")
            elif line.startswith("[ atoms ]"):
                read_atp = False
                read_itp = True
                read_atoms = True
            elif line.startswith("[ bonds ]"):
                read_atoms = False
            elif line.startswith("[ system ]"):
                read_atp = False
                read_itp = False
            if read_atp:
                atpf.write(line)
            if read_itp:
                itpf.write(line)
            if read_atoms and (not line.startswith(';')) and (not line.startswith('[ atoms ]')) and (line.strip()):
                atoms.append(line.split()[4])
    atpf.close()
    itpf.close()
    with open(posre, 'w') as f:
        f.write("[ position_restraints ]\n")
        f.write("; atom  type    fx    fy    fz\n")
        for i in range(len(atoms)):
            if not atoms[i].startswith("H"):
                f.write(f"{i+1:>6}     1  1000  1000 1000\n")


def run_command(cmd: List[str], **kwargs):
    cmd = [str(x) for x in cmd]
    sub = subprocess.Popen(
        args=cmd,
        stdin=subprocess.PIPE,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        **kwargs
    )
    out, err = sub.communicate()
    return_code = sub.poll()
    if sub.returncode != 0:
        raise ValueError(" ".join(cmd) + f" Failed to executed: \n {err.decode('utf-8')}")
    return return_code, out, err


def setup(wdir: os.PathLike, top: os.PathLike, gro: os.PathLike, config: os.PathLike):
    from gromacs.fileformats.mdp import MDP

    with open(config, 'r') as f:
        jdata = json.load(f)
    
    # prepare mdp settings
    mdp_settings = jdata['mdp']
    coul_lambdas = mdp_settings.get(
        "coul_lambdas", 
        [0.00, 0.25, 0.50, 0.75, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00]
    )
    vdw_lambdas = mdp_settings.get(
        "vdw_lambdas",
        [0.00, 0.00, 0.00, 0.00, 0.00, 0.05, 0.10, 0.20, 0.30, 0.40, 0.50, 0.60, 0.65, 0.70, 0.75, 0.80, 0.85, 0.90, 0.95, 1.00]
    )
    assert len(coul_lambdas) == len(vdw_lambdas)
    num_lambdas = len(coul_lambdas)
    stages = ["em", "nvt", "npt", "prod"]
    
    wdir = Path(wdir).resolve()
    wdir.mkdir(exist_ok=True, parents=True)

    # prepare topology & gro
    split_top(top, wdir / "MOL.atp", wdir / "MOL.itp", wdir / "posre_MOL.itp")
    shutil.copyfile(Path(__file__).parent / "template.top", wdir / "topol.top")
    shutil.copyfile(gro, wdir / "MOL.gro")

    if shutil.which("gmx_mpi") is not None:
        gmx = "gmx_mpi"
    elif shutil.which("gmx") is not None:
        gmx = "gmx"
    else:
        raise Exception("gmx not found")
    
    pwd = os.getcwd()
    os.chdir(wdir)
    run_command([gmx, "editconf", "-f", "MOL.gro", "-o", "newbox.gro", "-c", "-d", str(1.2), "-bt", "dodecahedron"])
    run_command([gmx, "solvate", "-cp", "newbox.gro", "-cs", "spc216.gro", "-o", "solv.gro", "-p", "topol.top"])
    os.chdir(pwd)

    for stage in stages:
        mdp = MDP(Path(__file__).parent / f"{stage}.mdp")
        mdp_settings[stage]['coul_lambdas'] = " ".join(f"{x:.2f}" for x in coul_lambdas)
        mdp_settings[stage]['vdw_lambdas'] = " ".join(f"{x:.2f}" for x in vdw_lambdas)
        for k, v in mdp_settings.get(stage, {}).items():
            mdp[k] = v
        for i in range(num_lambdas):
            mdp['init_lambda_state'] = i
            lambda_dir = wdir / f'lambda{i}'
            lambda_dir.mkdir(exist_ok=True, parents=True)
            script = lambda_dir / 'script.sh'
            if not script.is_file():
                shutil.copyfile(Path(__file__).parent / 'script.sh', script)
            stage_dir = lambda_dir / stage
            stage_dir.mkdir(exist_ok=True, parents=True)
            mdp.write(stage_dir / f'{stage}.mdp')
            if i == 0:
                run_command([
                    gmx, "grompp", 
                    "-c", str(wdir / "solv.gro"), 
                    "-f", str(stage_dir / f'{stage}.mdp'), 
                    '-p', str(wdir / 'topol.top'), 
                    '-r', str(wdir / 'solv.gro'),
                    '-pp', str(stage_dir / 'processed.top'), 
                    '-maxwarn', 10
                ])
            else:
                shutil.copyfile(wdir / f'lambda0/{stage}/processed.top', stage_dir / 'processed.top')    
            if stage == "em":
                shutil.copyfile(wdir / "solv.gro", stage_dir / "solv.gro")
    
    for file in Path.cwd().glob("*topol.tpr*"):
        file.unlink()

def run_md(wdir: os.PathLike, config: os.PathLike, exit_on_submit: bool = True):
    from dpdispatcher import Machine, Resources, Task, Submission

    with open(config, 'r') as f:
        jdata = json.load(f)
    
    machine = Machine.load_from_dict(jdata['machine'])
    resources = Resources.load_from_dict(jdata['resources'])
    num_lambdas = len(list(Path(wdir).glob("lambda*")))
    task_list = []
    for i in range(num_lambdas):
        task = Task(
            command=f"bash script.sh",
            task_work_path=f"lambda{i}",
            forward_files=[
                "script.sh",
                "em/processed.top", "nvt/processed.top", "npt/processed.top", "prod/processed.top", 
                "em/solv.gro", 
                "em/em.mdp", "nvt/nvt.mdp", "npt/npt.mdp", "prod/prod.mdp"
            ],
            backward_files=[
                "script.sh", "err", "log",
                "em/em.log", "nvt/nvt.log", "npt/npt.log", "prod/prod.log",
                "em/em.gro", "nvt/nvt.gro", "npt/npt.gro", "prod/prod.gro", "prod/prod.xtc", "prod/prod.xvg"
            ]
        )
        task_list.append(task)
    
    sub = Submission(work_base=str(wdir), machine=machine, resources=resources, task_list=task_list)
    res = sub.run_submission(exit_on_submit=exit_on_submit)
    return res


def analysis(wdir: os.PathLike):
    import alchemlyb
    from alchemlyb.parsing.gmx import extract_u_nk
    from alchemlyb.visualisation import plot_convergence, plot_mbar_overlap_matrix
    from alchemlyb.convergence import forward_backward_convergence
    from alchemlyb.estimators import MBAR

    wdir = Path(wdir).resolve()
    num_lambdas = len(list(wdir.glob("lambda*")))

    T = 298.15
    kBT = 8.314 * T / 1000 / 4.184 # kcal/mol

    u_nks = []
    for i in range(num_lambdas):
        xvg = str(wdir / f"lambda{i}/prod/prod.xvg")
        u_nks.append(extract_u_nk(xvg, T=T))
    # evaluate free energy with MBAR
    mbarEstimator = MBAR()
    mbarEstimator.fit(alchemlyb.concat(u_nks))
    dG = mbarEstimator.delta_f_.iloc[-1, 0] * kBT
    dG_std = mbarEstimator.d_delta_f_.iloc[-1, 0] * kBT
    # convergence analysis
    conv_df = forward_backward_convergence(u_nks, "mbar")
    for key in ['Forward', 'Forward_Error', 'Backward', 'Backward_Error']:
        conv_df[key] *= kBT
    conv_df.to_csv(wdir / "convergence.csv", index=None)
    conv_ax = plot_convergence(conv_df)
    conv_ax.set_ylabel("$\Delta G$ (kcal/mol)")
    conv_ax.set_title(f"Convergence Analysis")
    conv_ax.figure.savefig(str(wdir / "convergence.png"), dpi=300)
    # overlap matrix
    overlap_ax = plot_mbar_overlap_matrix(mbarEstimator.overlap_matrix)
    overlap_ax.figure.savefig(str(wdir / "overlap.png"), dpi=300)

    res = {"dG": dG, "dG_std": dG_std}
    with open(wdir / 'result.json', 'w') as f:
        json.dump(res, f)
    return dG, dG_std


if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description="Program to calculate hydration free energy")
    subparsers = parser.add_subparsers(title="Valid subcommands", dest="command")

    parser.add_argument(
        "-c", "--coordinate", 
        help="coordinate file, .gro format",
        dest="gro", 
        type=str
    )
    parser.add_argument(
        "-p", "--topology",   
        help="topology file, .top format", 
        dest="top",
        type=str
    )
    parser.add_argument(
        "-f", "--config",
        dest="config",
        help="json file to specify configs",
        type=str
    )
    parser.add_argument(
        "-w", "--working_dir",
        dest="wdir",
        help="working directory",
        type=str
    )
    parser.add_argument(
        "-e", "--exit_on_submit",
        dest='exit_on_submit',
        action="store_true",
        help="exist after submission"
    )

    anal_parser = subparsers.add_parser("analysis")
    anal_parser.add_argument("-w", "--working_dir", dest='wdir')

    args = parser.parse_args()
    if args.command is None:
        setup(args.wdir, args.top, args.gro, args.config)
        res = run_md(args.wdir, args.config, args.exit_on_submit)
        with open(os.path.join(args.wdir, 'submit.json'), 'w') as f:
            json.dump(res, f)
    else:
        analysis(args.wdir)
