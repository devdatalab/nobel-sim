/* See Novosad, Asher, Farquharson and Iljazi (2025) (https://paulnovosad.com/pdf/nobel-prizes.pdf), Sections 4.2 and B.5 for context.

Note that the variable names use the term IQ, but the concept here is 'latent scientific ability', which is likely correlated with IQ
 but includes other factors, e.g. conscientiousness, grit, etc.. We use iq as shorthand, but we are not talking about iq.

*/

/**********************************************/
/* program simulate_parents: runs the simulation   */
/**********************************************/
cap prog drop simulate_parents
prog def simulate_parents, rclass

  syntax, population(real) num_laureates(real) science_entry_share(real) nobel_error_rate(real) failure_rate(real) heritability(real) assortative(real) income_talent_genetic_corr(real) income_ability_corr(real) [ymax(real 1) sample_size(real 1000000)]

  qui {
    
    /* calculate additional parameters */
    
    /* calculate the population share getting nobel prizes */
    local top_share = `num_laureates' / `population'

    /* adjust for failure rate. if 25% fail, we need 5/4ths the number of top talents to generate N winners */
    local adjusted_top_share = `top_share' / (1 - `failure_rate')
    
    /* adjust for share of people going into the sciences. If 20% enter the sciences, need 5x N winners */
    local adjusted_top_share = `adjusted_top_share' / `science_entry_share'
    
    /* adjust for nobel committee error rate --- for our purposes, equivalent to scientific failure */
    local adjusted_top_share = `adjusted_top_share' / (1 - `nobel_error_rate')
    
    /* how many SDs above mean is the marginal laureate? */
    local sd_threshold = invnormal(1 - `adjusted_top_share')

    local sd_threshold = invnormal(1 - `adjusted_top_share')
  
    noi di "---"
    noi di "`adjusted_top_share' of the population are winners, vs. `top_share' without filters."
    noi di "The marginal nobel laureate has an ability score that is " %5.2f (`sd_threshold') " SD above the mean."
  
    /* start the simulation */
    clear
    set obs `sample_size'
    
    /* Generate mother and father genetic component of talent */
    /* Set genetic component variance such that phenotype will be N(0, 1) */
    gen gv_iq_dad = rnormal(0, sqrt(`heritability'))

    /* Incorporate assortative mating such that mom phenotype is N(0, 1) but has `assortative' correlation with dad */
    gen gv_iq_mom = sqrt(`heritability') * (`assortative' * (gv_iq_dad / sqrt(`heritability')) + sqrt(1 - `assortative' ^ 2) * rnormal())
    
    /* generate environmental component of dad's talent */
    gen ev_dad = sqrt(1 - `heritability') * rnormal()
    
    /* generate environmental component of mom's talent, again with correlation `assortative' with dad's */
    gen ev_mom = sqrt(1 - `heritability') * (`assortative' * (ev_dad / sqrt(1 - `heritability')) + sqrt(1 - `assortative' ^ 2) * rnormal())
    
    /* combine environment and genetic scores to get parent phenotype. These are N(0, 1) with the assortative correlation. */
    gen iq_dad = gv_iq_dad + ev_dad
    gen iq_mom = gv_iq_mom + ev_mom
    gen iq_parent_mean = (iq_mom + iq_dad) / 2
    
    /* generate father's "income capability" using an income-talent correlation, such that "income capability" is also N(0, 1) */
    /* this is the genetic predisposition to factors that contribute to high income production */
    gen income_capability_dad = `income_talent_genetic_corr' * iq_dad + sqrt(1 - `income_talent_genetic_corr' ^ 2) * rnormal()

    /* father's income is a product of his ability, plus random shocks */
    gen income_dad = income_capability_dad * `income_ability_corr' + sqrt(1 - `income_ability_corr' ^ 2) * rnormal()
    
    /* generate the genetic value for the child (assume 1 child per mom/dad pair) */
    
    /* adjust Mendelian variance for assortative mating */
    local mv = 0.5 * (1 - `assortative') * `heritability'
    
    /* Child gets average of parent talent, plus a random environmental component */
    gen gv_child = (gv_iq_mom + gv_iq_dad) / 2 + rnormal(0, sqrt(`mv'))
    
    /* note the noise term here is required. Genetically it comes from
    segregation and recombination.  Intuitively, without adding this
    noise, if each child's genetic value were simply the mean of the
    parents' values, genetic variance in the population would be
    constantly decreasing, which we don't see.
    
    This term scaling will also give us the right heritability estimate.
    */
    
    /* generate the phenotypic component for children, by adding environmental noise.
    Note in real life high ability parents might provide better environments, but we're
    building a null hypothesis with zero correlated environmental effects.
    */
    gen iq_child = gv_child + rnormal(0, sqrt(1 - `heritability'))
    /* check the distribution of child iq and parent-child correlations */
    
    /* child talent is mean 0 SD 1 */
    // sum iq_child
    
    /* review correlation between parent talents and parent/child talent to confirm they match assumptions */
    /* reg iq_child iq_mom
    reg iq_child iq_dad
    reg iq_child iq_mom iq_dad
    corr iq_child iq_dad
    corr iq_child iq_parent_mean */
  
    /* characterize parent income distribution of laureates */
    gen winner = iq_child > `sd_threshold'
  
    /* keep only the winners */
    keep if winner
  
    count if winner
    local num_winners = r(N)
  
    count if normal(income_dad) < 0.8
    local low_dad80 = r(N)
    local low_dad80_share: di %5.1f (`low_dad80' / `num_winners') * 100
    count if normal(income_dad) < 0.9
    local low_dad90 = r(N)
    local low_dad90_share: di %5.1f (`low_dad90' / `num_winners') * 100
  
    /* create strings with key results */
    noi di "Laureates are among the `num_winners' / " (`population'/1000000) "M most talented."
    local s_low_dad80_share = "  `low_dad80_share'% of fathers in bottom 80% (vs. empirical 15%)"
    local s_low_dad90_share = "  `low_dad90_share'% of fathers in bottom 90% (vs. empirical 28%)"
    return scalar bottom80 = `low_dad80_share'
    return scalar bottom90 = `low_dad90_share'
  
    /* manually transform father income into rank percentiles */
    /* can't rank in the data because we already truncated it. so do it analytically with the z-score */  
    gen income_rank_dad = normal(income_dad) * 100
  
    /* display information on key stats from this simulation */
    noi di "  `s_low_dad80_share'"
    noi di "  `s_low_dad90_share'"
  
    /* adjust upper bound and midpoints for histogram bins */
    replace income_rank_dad = 99.999 if income_rank_dad == 100
    egen graph_rank = cut(income_rank_dad), at(0(5)100)
    replace graph_rank = graph_rank + 2.5
    
    /* draw the histogram */
    local ymax = `ymax' * 100
    histogram graph_rank, discrete width(5) start(0) percent xtitle("Percentile") ///
      ytitle("Nobel Laureates (%)") xlabel(0(5)100) ylabel(0(5)`ymax') ytick(0(5)5) ///
      color(blue) lcolor(blue) lwidth(vthin) ///
      xtick(0(5)100) mcolor(navy) msize(medlarge) msymbol(o) ///
     	graphregion(fcolor(white))
  }
end
/** END program simulate_parents ***********/

/* run each sim 10 times and record the results */
cap erase $out/nobel_sims.csv
append_to_file using $out/nobel_sims.csv, s("i,sim_name,bottom80,bottom90") format(string) erase

global sample 1e8
forval i = 1/10 {
  di "-----------------"
  di "Iteration `i'/10:"

  /* nature-leaning scenario */
  simulate_parents, sample_size($sample) population(70000000) num_laureates(196) science_entry_share(0.1) nobel_error_rate(0.25) failure_rate(0.5) heritability(0.7) assortative(0.4) income_talent_genetic_corr(.9) income_ability_corr(0.6) ymax(.7)
  append_to_file using $out/nobel_sims.csv, s("`i',nature,`r(bottom80)',`r(bottom90)'") format(string) 
  graphout sim_rank_dad_nature, pdf
  
  /* mainline scenario */
  simulate_parents, sample_size($sample) population(70000000) num_laureates(196) science_entry_share(0.05) nobel_error_rate(0.35) failure_rate(0.75) heritability(0.5) assortative(0.3) income_talent_genetic_corr(.8) income_ability_corr(0.5) ymax(.7)
  append_to_file using $out/nobel_sims.csv, s("`i',mainline,`r(bottom80)',`r(bottom90)'") format(string) 
  graphout sim_rank_dad_mainline, pdf

  /* nurture-leaning */
  simulate_parents, sample_size($sample) population(70000000) num_laureates(196) science_entry_share(0.025) nobel_error_rate(0.45) failure_rate(0.85) heritability(0.3) assortative(0.2) income_talent_genetic_corr(.7) income_ability_corr(0.4) ymax(.7)
  append_to_file using $out/nobel_sims.csv, s("`i',nurture,`r(bottom80)',`r(bottom90)'") format(string) 
  graphout sim_rank_dad_mainline, pdf

  /* AI-parameter simulations */

  /* Gemini */
  simulate_parents, sample_size($sample) population(70000000) num_laureates(196) science_entry_share(0.2) nobel_error_rate(0.33) failure_rate(0.9) heritability(0.65) assortative(0.5) income_talent_genetic_corr(.6) income_ability_corr(0.5) ymax(.7)
  append_to_file using $out/nobel_sims.csv, s("`i',gemini,`r(bottom80)',`r(bottom90)'") format(string) 
  
  /* Claude */
  simulate_parents, sample_size($sample) population(70000000) num_laureates(196) science_entry_share(0.15) nobel_error_rate(0.25) failure_rate(0.85) heritability(0.6) assortative(0.4) income_talent_genetic_corr(.65) income_ability_corr(0.45) ymax(.7)
  append_to_file using $out/nobel_sims.csv, s("`i',claude,`r(bottom80)',`r(bottom90)'") format(string) 
  
  /* o3 */
  simulate_parents, sample_size($sample) population(70000000) num_laureates(196) science_entry_share(0.15) nobel_error_rate(0.35) failure_rate(0.97) heritability(0.55) assortative(0.30) income_talent_genetic_corr(.55) income_ability_corr(0.45) ymax(.7)
  append_to_file using $out/nobel_sims.csv, s("`i',o3,`r(bottom80)',`r(bottom90)'") format(string) 

  /* grok */
  simulate_parents, sample_size($sample) population(70000000) num_laureates(196) science_entry_share(0.5) nobel_error_rate(0.1) failure_rate(0.9) heritability(0.7) assortative(0.4) income_talent_genetic_corr(.67) income_ability_corr(0.45) ymax(.7)
  append_to_file using $out/nobel_sims.csv, s("`i',grok,`r(bottom80)',`r(bottom90)'") format(string) 
}

/* ---------------------------- cell: create simulation results table in latex-------------------- */

/* review simulation results */
import delimited using $out/nobel_sims.csv, clear varnames(1)

/* show simulation results on screen -- compare with 15% and 28% empirical */
tabstat bottom80 bottom90, by(sim_name) s(mean sd)

/* create simulation results stata-tex input file  */
global f $out/nobel_sim_results.csv
cap erase $f

/* loop over each simulation */
foreach sim in nature mainline nurture claude gemini o3 grok {
  foreach pct in 80 90 {

    if `pct' == 80 local empirical_pct 20.7
    if `pct' == 90 local empirical_pct 32.9
    
    /* summarize the bottom80 and bottom90 score */
    qui sum bottom`pct' if sim_name == "`sim'"

    /* store the mean and SD in the table_from_tpl input CSV */
    local mean: di %5.1f `r(mean)'
    local sd: di %5.1f `r(sd)'
    append_to_file using $f, s("`sim'_mean`pct',`mean'")
    append_to_file using $f, s("`sim'_sd`pct',`sd'")

    /* store share of explained inequality */
    local perc_explained: di %2.0f (100 * (`pct' - `mean') / (`pct' - `empirical_pct'))
    append_to_file using $f, s("`sim'_expl`pct',`perc_explained'")

    /* Write to the screen explaining how much was explained by genetics */
    di %10s "`sim'" ": `mean'% of laureates from bottom `pct': " %2.0f `perc_explained' "% of inequality was explained"
  
   }
}

/* generate the latex file with stata-tex */
table_from_tpl, t($pcode/tpl/sim_results_tpl.tex) r($f) o($out/sim_results.tex)

/* confirm it looks good */
cat $out/sim_results.tex


