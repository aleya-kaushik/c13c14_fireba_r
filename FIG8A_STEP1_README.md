# Figure 8A Step 1: restart-authoritative dynamic vegetation state

This branch is a deliberately narrow first change for the Figure 8A dynamic-PFT work.

Modified files:

- `restart_write.F90`
- `restart_read.F90`

Change:

- `lu_pftref` and `lu_area` are written from current runtime `sib%...%ipft/larea`.
- for `ic_mode='restart'`, those variables are restored into runtime state before `sibtype_setup`.
- fresh/requib initialization remains based on the external `sib_vs`.
- `subpref/sublarea` are synchronized when a restart is read so newly opened output files begin with metadata consistent with the restored vegetation state.

Scope limitation:

- `g_nlu` and LU allocation are still determined by the initial `sib_vs`; this supports an in-place one-LU PFT change but not yet creation of extra LUs for mixed PFTs.

No `tools_post` change is required for these two existing 2-D restart variables.

Before use in production, rebuild cleanly on the target HPC and rerun the established exact same-PFT restart regression.
