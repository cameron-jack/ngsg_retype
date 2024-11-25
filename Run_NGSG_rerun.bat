REM Don't auto update 
REM cmd /k "..\NGSgeno\venv_ngsgeno\Scripts\activate & git pull & streamlit run --server.enableXsrfProtection false ngsg_rerun.py"
cmd /k "..\NGSgeno\venv_ngsgeno\Scripts\activate & streamlit run --server.enableXsrfProtection false ngsg_rerun.py"