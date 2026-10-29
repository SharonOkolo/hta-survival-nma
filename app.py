"""Streamlit app: survival extrapolation and network meta-analysis results (computed in R, see analysis/)."""
import os
import pandas as pd
import plotly.express as px
import streamlit as st

OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "outputs")
r = lambda f: pd.read_csv(os.path.join(OUT, f))
st.set_page_config(page_title="Survival & NMA for HTA", layout="wide")

st.title("Survival extrapolation and network meta-analysis")
st.caption("Two analyses a health technology assessment needs, run in R (survival, flexsurv, netmeta, gemtc/JAGS). "
           "Data: colon cancer adjuvant trial (survival::colon) and 24 smoking-cessation trials (netmeta::smokingcessation).")

km, cox, ms, curves, kmp = r("km_summary.csv"), r("cox.csv"), r("mean_survival.csv"), r("extrapolation_curves.csv"), r("km_points.csv")
hr = cox[(cox.model == "unadjusted") & (cox.term == "rxLev+5FU")].iloc[0]
bay, freq, dic, sd = r("nma_bayesian.csv"), r("nma_frequentist.csv"), r("nma_bayes_model_fit.csv"), r("nma_bayes_sd.csv")

c1, c2, c3, c4 = st.columns(4)
c1.metric("HR death, Lev+5FU vs Obs", f"{hr.hr:.2f}", help=f"95% CI {hr.lcl:.2f} to {hr.ucl:.2f}")
c2.metric("5-year survival gain", f"{(km.surv_5y[km.arm == 'Lev+5FU'].iloc[0] - km.surv_5y[km.arm == 'Obs'].iloc[0]) * 100:.1f} pts")
c3.metric("Life-years gained (range)", f"{ms.life_years_gained.min():.1f} to {ms.life_years_gained.max():.1f}")
best = bay.loc[bay.sucra.idxmax()]
c4.metric("Best smoking intervention", best.treatment, f"SUCRA {best.sucra:.0%}", delta_color="off")

t1, t2, t3, t4 = st.tabs(["Survival extrapolation", "Subgroups", "Network meta-analysis", "Methods & validation"])

with t1:
    model = st.selectbox("Parametric model", ms.model.tolist(), index=ms.model.tolist().index("lognormal"))
    cv = curves[curves.model == model]
    fig = px.line(cv, x="time", y="surv", color="arm", labels={"time": "Years", "surv": "Overall survival"},
                  color_discrete_map={"Obs": "#9AA5B1", "Lev+5FU": "#2563EB"})
    for arm, colour in [("Obs", "#52606D"), ("Lev+5FU", "#1E3A8A")]:
        k = kmp[kmp.arm == arm]
        fig.add_scatter(x=k.time, y=k.surv, mode="lines", line=dict(shape="hv", color=colour, width=1.5), name=f"KM {arm}")
    fig.update_yaxes(tickformat=".0%", range=[0, 1])
    st.plotly_chart(fig, use_container_width=True)
    row = ms[ms.model == model].iloc[0]
    a, b, c = st.columns(3)
    a.metric("Mean survival, Obs (30 y)", f"{row.obs:.2f} years")
    b.metric("Mean survival, Lev+5FU (30 y)", f"{row.lev5fu:.2f} years")
    c.metric("Life-years gained", f"{row.life_years_gained:.2f}")
    st.markdown("The fitted curves all match the trial data well in the first 8 years, yet they disagree by "
                f"**{ms.life_years_gained.max() - ms.life_years_gained.min():.1f} life-years** in the long run. "
                "Choosing the model is a judgement on goodness of fit (AIC/BIC), visual fit and clinical plausibility, "
                "and the choice should be tested in scenario analyses.")
    st.dataframe(r("parametric_fit.csv").round(1), hide_index=True, use_container_width=True)

with t2:
    sg = r("subgroups.csv")
    sg["label"] = sg.subgroup + ": " + sg.level
    fig = px.scatter(sg, x="hr", y="label", error_x=sg.ucl - sg.hr, error_x_minus=sg.hr - sg.lcl, log_x=True,
                     labels={"hr": "Hazard ratio, Lev+5FU vs Obs", "label": ""})
    fig.add_vline(x=1, line_dash="dash", line_color="grey")
    st.plotly_chart(fig, use_container_width=True)
    st.markdown("Exploratory and not adjusted for multiplicity. The sex interaction (p = "
                f"{sg.p_interaction[sg.subgroup == 'Sex'].iloc[0]:.2f}) is hypothesis-generating: on its own it shouldn't be used to restrict treatment.")
    st.dataframe(sg.drop(columns="label").round(3), hide_index=True, use_container_width=True)

with t3:
    comp = r("nma_compare.csv")
    fig = px.scatter(comp, x="or", y="treatment", color="method", error_x=comp.ucl - comp["or"], error_x_minus=comp["or"] - comp.lcl,
                     log_x=True, labels={"or": "Odds ratio of quitting vs no contact", "treatment": ""},
                     color_discrete_map={"Bayesian (gemtc/JAGS)": "#2563EB", "Frequentist (netmeta)": "#9AA5B1"})
    fig.add_vline(x=1, line_dash="dash", line_color="grey")
    st.plotly_chart(fig, use_container_width=True)
    c1, c2 = st.columns(2)
    c1.markdown("**Rank probabilities (Bayesian)**")
    c1.dataframe(r("nma_rank_probabilities.csv").round(2), hide_index=True, use_container_width=True)
    c2.markdown("**Model fit**")
    c2.dataframe(dic.round(3), hide_index=True, use_container_width=True)
    st.markdown("**League table (frequentist, random effects; row vs column OR)**")
    st.dataframe(r("nma_league_table.csv"), hide_index=True, use_container_width=True)
    st.markdown("**Inconsistency: node-splitting, direct vs indirect**")
    st.dataframe(r("nma_nodesplit.csv").round(3), hide_index=True, use_container_width=True)

with t4:
    st.markdown(open(os.path.join(os.path.dirname(os.path.abspath(__file__)), "SAP.md")).read())
