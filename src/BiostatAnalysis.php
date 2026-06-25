<?php

declare(strict_types=1);

namespace TouilElhadj\BiostatPhp;

use InvalidArgumentException;

/**
 * BiostatAnalysis — pure-PHP biostatistics library.
 *
 * Implements descriptive, bivariate and multivariate methods of inference
 * for survey-based epidemiological studies, with a feature set inspired by
 * the R packages `stats`, `car`, `lme4`, `geepack` and `mice`.
 *
 * **Method families**
 *
 *  • Descriptive:        {@see mean}, {@see std}, {@see median}, {@see quantile}
 *  • 2 × 2 tables:       {@see chi2Test2x2}, {@see oddsRatio}
 *  • Means:              {@see tTest} (Welch), {@see anova}
 *  • Correlation:        {@see pearson}, {@see spearman}
 *  • Binomial:           {@see binomialTest}
 *  • Logistic:           {@see logisticRegression}, {@see logisticRegressionMulti},
 *                        {@see hosmerLemeshow}
 *  • Multicollinearity:  {@see vif}
 *  • Logit linearity:    {@see boxTidwell}
 *  • Mixed models:       {@see glmmLogistic} (PQL),
 *                        {@see geeLogistic}  (GEE + Liang–Zeger sandwich)
 *  • Missing data:       {@see mice}, {@see rubinPool}
 *
 * **Verification**: every public method is cross-checked against R 4.x and
 * IBM SPSS Statistics 25; see `docs/validation-tables.md` for the numerical
 * comparison and `tests/` for the executable assertions.
 *
 * @author  Elhadj TOUIL <touilelhadj@live.com>
 * @license MIT
 */
class BiostatAnalysis
{
    use LinearAlgebra;
    use Distributions;

    /**
     * Optional row-oriented dataset that descriptive methods can pull
     * columns from with array_column().
     *
     * @var array<int, array<string, mixed>>
     */
    private array $data;

    /**
     * Default α used by inferential helpers when none is supplied.
     */
    private float $alpha = 0.05;

    /**
     * @param array<int, array<string, mixed>> $data row-oriented dataset.
     *        Each row is an associative array (e.g. as returned by
     *        PDOStatement::fetchAll(PDO::FETCH_ASSOC)). May be empty if
     *        only the raw numeric routines are needed.
     */
    public function __construct(array $data = [])
    {
        $this->data = $data;
    }

    // ═══════════════════════════════════════════════════════════════
    // DESCRIPTIVE STATISTICS
    // ═══════════════════════════════════════════════════════════════


    /**
     * Moyenne arithmétique
     */
    public function mean($array) {
        $array = array_filter($array, fn($v) => $v !== null && $v !== '');
        return count($array) > 0 ? array_sum($array) / count($array) : 0;
    }
    
    /**
     * Écart-type (n-1)
     */
    public function std($array) {
        $array = array_values(array_filter($array, fn($v) => is_numeric($v)));
        if(count($array) < 2) return 0;
        
        $mean = $this->mean($array);
        $variance = array_sum(array_map(fn($v) => ($v - $mean) ** 2, $array)) / (count($array) - 1);
        return sqrt($variance);
    }
    
    /**
     * Médiane
     */
    public function median($array) {
        $array = array_values(array_filter($array, fn($v) => is_numeric($v)));
        if(empty($array)) return 0;
        
        sort($array);
        $n = count($array);
        return $n % 2 ? $array[intdiv($n, 2)] : ($array[$n/2 - 1] + $array[$n/2]) / 2;
    }
    
    /**
     * Quantiles
     */
    public function quantile($array, $q) {
        $array = array_values(array_filter($array, fn($v) => is_numeric($v)));
        if(empty($array)) return 0;
        
        sort($array);
        $n = count($array);
        $index = ($n - 1) * $q;
        $lower = floor($index);
        $upper = ceil($index);
        
        if($lower == $upper) return $array[$lower];
        return $array[$lower] * ($upper - $index) + $array[$upper] * ($index - $lower);
    }
    
    // ═══════════════════════════════════════════════════════════════
    // TESTS STATISTIQUES INFÉRENTIELS
    // ═══════════════════════════════════════════════════════════════
    
    /**
     * Test Chi² d'indépendance (tableau 2×2) avec correction de continuité Yates par défaut
     * Aligné sur R chisq.test(..., correct=TRUE).
     * Retourne: ['chi2' => float, 'p' => float, 'df' => int]
     */
    public function chi2Test2x2($a, $b, $c, $d, $yates = true) {
        $n = $a + $b + $c + $d;
        if($n == 0 || ($a+$b) == 0 || ($c+$d) == 0 || ($a+$c) == 0 || ($b+$d) == 0) {
            return ['chi2' => 0, 'p' => 1, 'df' => 1];
        }
        // Yates : |ad-bc| - n/2 (clamped à 0)
        $diff = abs($a*$d - $b*$c);
        if($yates) $diff = max(0.0, $diff - $n/2.0);
        $chi2 = $n * pow($diff, 2) / (($a+$b) * ($c+$d) * ($a+$c) * ($b+$d));
        $p = $this->chi2CDF($chi2, 1);
        return [
            'chi2' => round($chi2, 4),
            'p' => round($p, 4),
            'df' => 1,
            'significant' => $p < $this->alpha
        ];
    }
    
    /**
     * Odds Ratio avec IC95%
     */
    /**
     * Odds ratio with 95 % Wald confidence interval.
     *
     * If any cell is zero, the Haldane–Anscombe continuity correction
     * (+0.5 added to every cell) is automatically applied to keep both
     * the point estimate and the CI defined.
     *
     * Reference for the correction:
     *   Anscombe, F. J. (1956). On estimating binomial response relations.
     *   Biometrika 43, 461–464.
     *
     * @param int|float $a exposed cases
     * @param int|float $b exposed non-cases
     * @param int|float $c unexposed cases
     * @param int|float $d unexposed non-cases
     *
     * @return array{
     *     or: float,
     *     ci_low: float|null,
     *     ci_high: float|null,
     *     correction_applied: bool
     * }
     */
    public function oddsRatio($a, $b, $c, $d) {
        // Haldane–Anscombe continuity correction when any cell is zero
        $correction = ($a == 0 || $b == 0 || $c == 0 || $d == 0);
        if ($correction) {
            $a += 0.5;
            $b += 0.5;
            $c += 0.5;
            $d += 0.5;
        }

        $or = ($a * $d) / ($b * $c);
        $se = sqrt(1 / $a + 1 / $b + 1 / $c + 1 / $d);
        $ci_low  = exp(log($or) - 1.96 * $se);
        $ci_high = exp(log($or) + 1.96 * $se);

        return [
            'or'                 => round($or, 2),
            'ci_low'             => round($ci_low,  2),
            'ci_high'            => round($ci_high, 2),
            'correction_applied' => $correction,
        ];
    }
    
    /**
     * Test t de Student (Welch)
     * Retourne: ['t' => float, 'p' => float, 'df' => float, 'm1' => float, 'm2' => float]
     */
    public function tTest($group1, $group2) {
        $g1 = array_values(array_filter($group1, 'is_numeric'));
        $g2 = array_values(array_filter($group2, 'is_numeric'));
        
        $n1 = count($g1);
        $n2 = count($g2);
        
        if($n1 < 2 || $n2 < 2) {
            return ['t' => 0, 'p' => 1, 'df' => 0, 'm1' => 0, 'm2' => 0];
        }
        
        $m1 = $this->mean($g1);
        $m2 = $this->mean($g2);
        $s1 = $this->std($g1);
        $s2 = $this->std($g2);
        
        // Erreur standard
        $se = sqrt(($s1**2)/$n1 + ($s2**2)/$n2);
        
        if($se == 0) {
            return ['t' => 0, 'p' => 1, 'df' => 0, 'm1' => $m1, 'm2' => $m2];
        }
        
        // Statistique t
        $t = ($m1 - $m2) / $se;
        
        // Degrés de liberté Welch
        $df = pow(($s1**2/$n1 + $s2**2/$n2), 2) / 
              (pow($s1**2/$n1, 2)/max($n1-1, 1) + pow($s2**2/$n2, 2)/max($n2-1, 1));
        
        // p-value approximative (bilatéral)
        $p = 2 * $this->studentTCDF(abs($t), $df);
        
        return [
            't' => round($t, 3),
            'p' => round($p, 4),
            'df' => round($df, 1),
            'm1' => round($m1, 2),
            'm2' => round($m2, 2),
            'sd1' => round($s1, 2),
            'sd2' => round($s2, 2),
            'n1' => $n1,
            'n2' => $n2,
            'significant' => $p < $this->alpha
        ];
    }
    
    /**
     * Corrélation de Pearson
     */
    public function pearson($x, $y) {
        // Filtrer paires valides
        $pairs = [];
        foreach($x as $i => $xi) {
            if(isset($y[$i]) && is_numeric($xi) && is_numeric($y[$i])) {
                $pairs[] = [(float)$xi, (float)$y[$i]];
            }
        }
        
        $n = count($pairs);
        if($n < 3) return ['r' => 0, 'p' => 1, 'n' => $n];
        
        // Extraire x et y
        $xs = array_column($pairs, 0);
        $ys = array_column($pairs, 1);
        
        $mx = $this->mean($xs);
        $my = $this->mean($ys);
        
        // Covariance et écarts-types
        $cov = 0;
        $ssx = 0;
        $ssy = 0;
        
        foreach($pairs as $p) {
            $dx = $p[0] - $mx;
            $dy = $p[1] - $my;
            $cov += $dx * $dy;
            $ssx += $dx * $dx;
            $ssy += $dy * $dy;
        }
        
        $r = ($ssx * $ssy > 0) ? $cov / sqrt($ssx * $ssy) : 0;
        
        // Test de significativité
        $t = abs($r) * sqrt($n - 2) / sqrt(max(1 - $r**2, 0.0001));
        $p = 2 * $this->studentTCDF($t, $n - 2);
        
        return [
            'r' => round($r, 3),
            'p' => round($p, 4),
            'n' => $n,
            'r2' => round($r**2, 3),
            'significant' => $p < $this->alpha
        ];
    }
    
    /**
     * Corrélation de Spearman (rangs)
     */
    public function spearman($x, $y) {
        // Convertir en rangs
        $rx = $this->ranks($x);
        $ry = $this->ranks($y);
        
        // Appliquer Pearson sur les rangs
        return $this->pearson($rx, $ry);
    }
    
    /**
     * ANOVA à un facteur
     */
    public function anova($groups) {
        // Aplatir tous les groupes
        $all = [];
        $ns = [];
        $means = [];
        
        foreach($groups as $key => $group) {
            $g = array_filter($group, 'is_numeric');
            $all = array_merge($all, $g);
            $ns[$key] = count($g);
            $means[$key] = $this->mean($g);
        }
        
        $n_total = count($all);
        $k = count($groups);
        
        if($n_total < $k + 1) {
            return ['F' => 0, 'p' => 1, 'dfB' => 0, 'dfW' => 0];
        }
        
        $grand_mean = $this->mean($all);
        
        // Somme des carrés inter-groupes (between)
        $ssB = 0;
        foreach($groups as $key => $group) {
            $ssB += $ns[$key] * pow($means[$key] - $grand_mean, 2);
        }
        
        // Somme des carrés intra-groupes (within)
        $ssW = 0;
        foreach($groups as $group) {
            $g = array_filter($group, 'is_numeric');
            $m = $this->mean($g);
            foreach($g as $v) {
                $ssW += pow($v - $m, 2);
            }
        }
        
        $dfB = $k - 1;
        $dfW = $n_total - $k;
        
        if($dfW == 0 || $ssW == 0) {
            return ['F' => 0, 'p' => 1, 'dfB' => $dfB, 'dfW' => $dfW];
        }
        
        $msB = $ssB / $dfB;
        $msW = $ssW / $dfW;
        
        $F = $msB / $msW;
        
        // p-value approximative
        $p = $this->fCDF($F, $dfB, $dfW);
        
        return [
            'F' => round($F, 3),
            'p' => round($p, 4),
            'dfB' => $dfB,
            'dfW' => $dfW,
            'means' => $means,
            'n' => $ns,
            'significant' => $p < $this->alpha
        ];
    }
    
    /**
     * Test binomial exact (proportion)
     */
    public function binomialTest($successes, $n, $p0 = 0.5) {
        if($n == 0) return ['p' => 1, 'ci_low' => 0, 'ci_high' => 1];
        
        $p_obs = $successes / $n;
        
        // IC95% Wilson
        $z = 1.96;
        $denom = 1 + $z**2 / $n;
        $center = ($p_obs + $z**2 / (2*$n)) / $denom;
        $margin = $z * sqrt($p_obs * (1 - $p_obs) / $n + $z**2 / (4*$n**2)) / $denom;
        
        $ci_low = max(0, $center - $margin);
        $ci_high = min(1, $center + $margin);
        
        // Test exact impossible en PHP pur, utiliser approximation normale
        $se = sqrt($p0 * (1 - $p0) / $n);
        $z_stat = abs($p_obs - $p0) / $se;
        $p = 2 * (1 - $this->normalCDF($z_stat));
        
        return [
            'p_obs' => round($p_obs, 4),
            'p' => round($p, 4),
            'ci_low' => round($ci_low, 4),
            'ci_high' => round($ci_high, 4),
            'significant' => $p < $this->alpha
        ];
    }
    
    /**
     * Régression logistique simple
     * Utilise la méthode de Newton-Raphson
     */
    public function logisticRegression($y, $x) {
        // Préparer données
        $n = count($y);
        $y = array_values($y);
        $x = array_values($x);
        
        // Initialisation
        $beta = [0, 0]; // [intercept, slope]
        $max_iter = 100;
        $tolerance = 0.0001;
        
        for($iter = 0; $iter < $max_iter; $iter++) {
            // Calculer prédictions et gradient
            $gradient = [0, 0];
            $hessian = [[0, 0], [0, 0]];
            
            for($i = 0; $i < $n; $i++) {
                $eta = $beta[0] + $beta[1] * $x[$i];
                $p = 1 / (1 + exp(-$eta));
                $w = $p * (1 - $p);
                
                $residual = $y[$i] - $p;
                
                $gradient[0] += $residual;
                $gradient[1] += $residual * $x[$i];
                
                $hessian[0][0] -= $w;
                $hessian[0][1] -= $w * $x[$i];
                $hessian[1][0] -= $w * $x[$i];
                $hessian[1][1] -= $w * $x[$i] * $x[$i];
            }
            
            // Inverser Hessian
            $det = $hessian[0][0] * $hessian[1][1] - $hessian[0][1] * $hessian[1][0];
            if(abs($det) < 1e-10) break;
            
            $inv_hess = [
                [$hessian[1][1] / $det, -$hessian[0][1] / $det],
                [-$hessian[1][0] / $det, $hessian[0][0] / $det]
            ];
            
            // Mise à jour
            $delta = [
                $inv_hess[0][0] * $gradient[0] + $inv_hess[0][1] * $gradient[1],
                $inv_hess[1][0] * $gradient[0] + $inv_hess[1][1] * $gradient[1]
            ];
            
            $beta[0] -= $delta[0];
            $beta[1] -= $delta[1];
            
            // Convergence?
            if(abs($delta[0]) < $tolerance && abs($delta[1]) < $tolerance) break;
        }
        
        // Calculer OR et IC95%
        $or = exp($beta[1]);
        $se = sqrt(-1 / $hessian[1][1]); // Erreur standard
        $ci_low = exp($beta[1] - 1.96 * $se);
        $ci_high = exp($beta[1] + 1.96 * $se);
        
        // p-value (Wald)
        $z = $beta[1] / $se;
        $p = 2 * (1 - $this->normalCDF(abs($z)));
        
        return [
            'intercept' => round($beta[0], 4),
            'coef' => round($beta[1], 4),
            'or' => round($or, 2),
            'ci_low' => round($ci_low, 2),
            'ci_high' => round($ci_high, 2),
            'p' => round($p, 4),
            'significant' => $p < $this->alpha
        ];
    }
    
    // ═══════════════════════════════════════════════════════════════
    // FONCTIONS UTILITAIRES
    // ═══════════════════════════════════════════════════════════════
    
    /**
     * Convertir en rangs (pour Spearman)
     */
    private function ranks($array) {
        $array = array_filter($array, fn($v) => is_numeric($v));
        $sorted = $array;
        asort($sorted);
        
        $ranks = [];
        $rank = 1;
        foreach($sorted as $key => $value) {
            $ranks[$key] = $rank++;
        }
        
        return $ranks;
    }
    
    
    // ═══════════════════════════════════════════════════════════════
    // RÉGRESSION LOGISTIQUE MULTIVARIÉE (Newton-Raphson généralisé)
    // ═══════════════════════════════════════════════════════════════

    /**
     * Régression logistique multivariée
     * @param array $y      Variable dépendante binaire (0/1), n obs
     * @param array $X      Matrice design n × k (sans intercept — ajouté auto)
     * @param array $names  Noms des k variables (pour le rapport)
     * @return array        ['coef'=>[], 'se'=>[], 'or'=>[], 'ci_low'=>[], 'ci_high'=>[],
     *                       'p'=>[], 'aic'=>x, 'auc'=>x, 'hl'=>['chi2'=>x,'p'=>x],
     *                       'converged'=>bool, 'iter'=>n, 'predicted'=>[]]
     */
    public function logisticRegressionMulti($y, $X, $names = []) {
        $n = count($y);
        $k = count($X[0] ?? []);
        if($n < 10 || $k < 1) return ['error'=>'insufficient data'];

        // Ajouter colonne intercept (1) en première position
        $design = [];
        for($i=0; $i<$n; $i++){
            $row = [1.0];
            for($j=0; $j<$k; $j++) $row[] = (float)$X[$i][$j];
            $design[] = $row;
        }
        $K = $k + 1;
        $beta = array_fill(0, $K, 0.0);
        $max_iter = 100;
        $tol = 1e-6;
        $converged = false; $iter = 0;

        for($it=0; $it<$max_iter; $it++){
            $iter = $it+1;
            // Compute predictions, gradient and Hessian
            $grad = array_fill(0, $K, 0.0);
            $hess = [];
            for($a=0;$a<$K;$a++) $hess[$a]=array_fill(0,$K,0.0);

            for($i=0;$i<$n;$i++){
                $eta = 0.0;
                for($j=0;$j<$K;$j++) $eta += $beta[$j] * $design[$i][$j];
                $p = 1.0 / (1.0 + exp(-max(min($eta,30),-30)));
                $w = $p*(1-$p);
                $resid = ((float)$y[$i]) - $p;
                for($j=0;$j<$K;$j++){
                    $grad[$j] += $resid * $design[$i][$j];
                    for($l=0;$l<$K;$l++) $hess[$j][$l] -= $w * $design[$i][$j] * $design[$i][$l];
                }
            }

            // Inverse Hessian (négative définie) — invertir matrice K×K
            $inv = $this->matrixInverse($hess);
            if($inv === null) break;

            // Update : beta -= inv(H) · grad
            $delta = array_fill(0,$K,0.0);
            for($a=0;$a<$K;$a++)
                for($b=0;$b<$K;$b++)
                    $delta[$a] += $inv[$a][$b] * $grad[$b];

            $maxd = 0;
            for($a=0;$a<$K;$a++){ $beta[$a] -= $delta[$a]; if(abs($delta[$a])>$maxd) $maxd=abs($delta[$a]); }
            if($maxd < $tol){ $converged=true; break; }
        }

        // Erreurs standard à partir de inv(-Hessian)
        $invH = $this->matrixInverse($hess);
        $se = array_fill(0,$K,null);
        if($invH !== null){
            for($j=0;$j<$K;$j++){
                $v = -$invH[$j][$j];
                $se[$j] = $v>0 ? sqrt($v) : null;
            }
        }

        // OR + IC95 + p
        $coef = $beta;
        $or = []; $cilow=[]; $cihigh=[]; $pval=[];
        for($j=0;$j<$K;$j++){
            $or[$j] = exp($coef[$j]);
            if($se[$j] !== null){
                $cilow[$j]  = exp($coef[$j] - 1.96*$se[$j]);
                $cihigh[$j] = exp($coef[$j] + 1.96*$se[$j]);
                $z = $coef[$j] / $se[$j];
                $pval[$j] = 2 * (1 - $this->normalCDF(abs($z)));
            } else {
                $cilow[$j]=null; $cihigh[$j]=null; $pval[$j]=null;
            }
        }

        // Prédictions p_i
        $predicted = [];
        for($i=0;$i<$n;$i++){
            $eta=0; for($j=0;$j<$K;$j++) $eta += $beta[$j]*$design[$i][$j];
            $predicted[$i] = 1/(1+exp(-max(min($eta,30),-30)));
        }

        // Log-likelihood pour AIC
        $ll = 0;
        for($i=0;$i<$n;$i++){
            $p = max(min($predicted[$i], 0.9999), 0.0001);
            $ll += $y[$i]*log($p) + (1-$y[$i])*log(1-$p);
        }
        $aic = 2*$K - 2*$ll;

        // AUC (Hanley-McNeil via Mann-Whitney)
        $auc = $this->aucFromPredictions($y, $predicted);

        // Hosmer-Lemeshow goodness-of-fit
        $hl = $this->hosmerLemeshow($y, $predicted, 10);

        // Bâtir résultat lisible : exclure intercept dans tableaux secondaires
        $varNames = ['(Intercept)'];
        for($j=0;$j<$k;$j++) $varNames[] = $names[$j] ?? "X$j";

        $result = [
            'converged'=>$converged, 'iter'=>$iter, 'n'=>$n, 'k'=>$k,
            'names'    => $varNames,
            'coef'     => array_map(fn($v)=>round($v,4), $coef),
            'se'       => array_map(fn($v)=>$v===null?null:round($v,4), $se),
            'or'       => array_map(fn($v)=>round($v,3), $or),
            'ci_low'   => array_map(fn($v)=>$v===null?null:round($v,3), $cilow),
            'ci_high'  => array_map(fn($v)=>$v===null?null:round($v,3), $cihigh),
            'p'        => array_map(fn($v)=>$v===null?null:round($v,4), $pval),
            'aic'      => round($aic, 2),
            'log_lik'  => round($ll, 2),
            'auc'      => $auc===null?null:round($auc, 3),
            'hl'       => $hl,
            'predicted'=> array_map(fn($v)=>round($v,4), $predicted),
        ];
        return $result;
    }

    /* Inversion matrice carrée — Gauss-Jordan */

    /* AUC à partir des prédictions (Mann-Whitney) */
    private function aucFromPredictions($y, $p){
        $n=count($y);
        $pos=[]; $neg=[];
        for($i=0;$i<$n;$i++) ($y[$i]==1) ? $pos[]=$p[$i] : $neg[]=$p[$i];
        $np = count($pos); $nn = count($neg);
        if($np==0||$nn==0) return null;
        $sum = 0;
        foreach($pos as $a) foreach($neg as $b){
            if($a>$b) $sum++;
            elseif(abs($a-$b)<1e-9) $sum += 0.5;
        }
        return $sum / ($np*$nn);
    }

    /**
     * Hosmer-Lemeshow goodness-of-fit
     * H0 : le modèle est bien ajusté (p>0.05 souhaité).
     */
    public function hosmerLemeshow($y, $p, $g=10){
        $n = count($y);
        if($n<$g) return ['chi2'=>null,'p'=>null,'df'=>null];
        // Trier par p croissant
        $pairs=[];
        for($i=0;$i<$n;$i++) $pairs[] = ['p'=>$p[$i],'y'=>$y[$i]];
        usort($pairs, fn($a,$b)=>$a['p']<=>$b['p']);
        // Découper en g groupes de taille ≈ équivalente
        $size = intdiv($n,$g);
        $chi2=0; $valid_groups=0;
        for($k=0;$k<$g;$k++){
            $start=$k*$size;
            $end=($k===$g-1)? $n : ($k+1)*$size;
            $obs=0; $exp_=0;
            for($i=$start;$i<$end;$i++){
                $obs += $pairs[$i]['y'];
                $exp_ += $pairs[$i]['p'];
            }
            $m = $end - $start;
            if($m===0 || $exp_<=0 || $exp_>=$m) continue;
            $chi2 += pow($obs-$exp_,2)/($exp_*(1-$exp_/$m));
            $valid_groups++;
        }
        $df = max(1, $valid_groups - 2);
        $p_hl = $this->chi2CDF($chi2, $df);
        return ['chi2'=>round($chi2,3),'p'=>round($p_hl,4),'df'=>$df];
    }

    /**
     * Correction Benjamini-Hochberg pour tests multiples
     * @param array $pvalues p-values brutes (clés conservées)
     * @return array         p-values ajustées (FDR)
     */
    public static function benjaminiHochberg(array $pvalues): array {
        $valid = array_filter($pvalues, fn($p)=>is_numeric($p) && $p>=0 && $p<=1);
        if(empty($valid)) return $pvalues;
        $keys = array_keys($valid);
        $vals = array_values($valid);
        $m = count($vals);
        // Trier par p croissant
        $idx = range(0,$m-1);
        usort($idx, fn($a,$b)=>$vals[$a] <=> $vals[$b]);
        $adj = array_fill(0,$m,1.0);
        $prev = 1.0;
        for($i=$m-1;$i>=0;$i--){
            $rank = $i+1;
            $bh = $vals[$idx[$i]] * $m / $rank;
            $bh = min(1.0, $bh, $prev);
            $adj[$idx[$i]] = $bh;
            $prev = $bh;
        }
        $out = $pvalues;
        for($i=0;$i<$m;$i++) $out[$keys[$i]] = round($adj[$i],4);
        return $out;
    }

    /**
     * Hosmer-Lemeshow utilise la chi2CDF existante (privée) — on la rend
     * accessible en interne via $this->chi2CDF($chi2,$df)
     */

    // ═══════════════════════════════════════════════════════════════
    // VIF — VARIANCE INFLATION FACTOR
    // ═══════════════════════════════════════════════════════════════
    /**
     * Calcule le VIF pour chaque colonne prédictrice continue de X.
     * VIF_j = 1 / (1 - R²_j), où R²_j provient de la régression OLS
     * de X_j sur toutes les autres colonnes (intercept inclus).
     * Seuil de vigilance courant : VIF > 2.5 (Allison, 2012)
     *                              VIF > 5 ou 10 (seuils alternatifs)
     *
     * @param array $X       n×k matrice prédicteurs (sans intercept)
     * @param array $names   nom de chaque variable
     * @return array         ['VarName' => VIF_value, ...]
     */
    public function vif(array $X, array $names = []): array {
        $n = count($X);
        if($n < 10) return ['error' => 'insufficient data'];
        $k = count($X[0] ?? []);
        if($k < 2) return ['error' => 'VIF requires >= 2 predictors'];

        $out = [];
        for($j = 0; $j < $k; $j++){
            // Cible = colonne j
            $y_j = [];
            for($i = 0; $i < $n; $i++) $y_j[] = (float)$X[$i][$j];

            // Régresseurs = autres colonnes + intercept
            $X_others = [];
            for($i = 0; $i < $n; $i++){
                $row = [1.0]; // intercept
                for($l = 0; $l < $k; $l++) if($l !== $j) $row[] = (float)$X[$i][$l];
                $X_others[] = $row;
            }
            $reg = $this->olsRegression($y_j, $X_others);
            if($reg === null){
                $out[$names[$j] ?? "X$j"] = ['vif' => null, 'r2' => null, 'error' => 'singular'];
                continue;
            }
            $r2 = $reg['r2'];
            $vif_val = ($r2 < 0.9999) ? (1.0 / (1.0 - $r2)) : INF;
            $out[$names[$j] ?? "X$j"] = [
                'vif'   => is_finite($vif_val) ? round($vif_val, 3) : 'Inf',
                'r2'    => round($r2, 4),
                'flag'  => $vif_val > 5    ? 'CRITICAL'
                        : ($vif_val > 2.5 ? 'WARNING' : 'OK'),
                'tolerance' => round(1 - $r2, 4), // = 1/VIF
            ];
        }
        return $out;
    }

    // ═══════════════════════════════════════════════════════════════
    // BOX-TIDWELL — TEST DE LINÉARITÉ DU LOGIT
    // ═══════════════════════════════════════════════════════════════
    /**
     * Test de Box-Tidwell (1962) : pour chaque prédicteur continu X_j
     * (qui DOIT être strictement positif), ajoute le terme d'interaction
     * X_j × ln(X_j) au modèle logistique multivarié.
     * H0 : le logit est linéaire en X_j (β du terme d'interaction = 0).
     * Rejet de H0 (p<0.05) → non-linéarité ; recoder en classes ou splines.
     *
     * @param array $y                 outcome binaire
     * @param array $X                 n×k matrice prédicteurs (sans intercept)
     * @param array $continuous_idx    indices (0-based) des colonnes continues à tester
     * @param array $names             noms des prédicteurs
     * @return array                   [Var => ['p'=>, 'beta_xlnx'=>, 'verdict'=>]]
     */
    public function boxTidwell(array $y, array $X, array $continuous_idx, array $names = []): array {
        $results = [];
        $n = count($X);
        $k = count($X[0] ?? []);

        foreach($continuous_idx as $j){
            // Vérifier que toutes les valeurs sont > 0
            $shift = 0.0;
            $minv  = INF;
            for($i = 0; $i < $n; $i++){
                if($X[$i][$j] < $minv) $minv = $X[$i][$j];
            }
            if($minv <= 0){
                // Shift positif (Box-Tidwell exige X > 0)
                $shift = abs($minv) + 1.0;
            }
            // Construire X augmenté avec colonne X_j × ln(X_j + shift)
            $X_aug = [];
            for($i = 0; $i < $n; $i++){
                $row = [];
                for($l = 0; $l < $k; $l++) $row[] = (float)$X[$i][$l];
                $xj = (float)$X[$i][$j] + $shift;
                $row[] = $xj * log($xj);
                $X_aug[] = $row;
            }

            $name_j = $names[$j] ?? "X$j";
            $names_aug = $names;
            $names_aug[] = $name_j . "×ln";

            $fit = $this->logisticRegressionMulti($y, $X_aug, $names_aug);
            if(isset($fit['error']) || !$fit['converged']){
                $results[$name_j] = ['error' => 'fit failed', 'p' => null];
                continue;
            }

            $K = count($fit['coef']);
            $p_val = $fit['p'][$K - 1];
            $beta_int = $fit['coef'][$K - 1];
            $se_int = $fit['se'][$K - 1];

            $results[$name_j] = [
                'beta_xlnx'  => round($beta_int, 4),
                'se'         => $se_int === null ? null : round($se_int, 4),
                'p'          => $p_val === null ? null : round($p_val, 4),
                'shift'      => round($shift, 3),
                'linear_ok'  => ($p_val !== null && $p_val > 0.05),
                'verdict'    => ($p_val !== null && $p_val > 0.05)
                                ? 'Linéarité du logit non rejetée'
                                : ($p_val === null ? 'Indéterminé' : 'Non-linéarité détectée'),
            ];
        }
        return $results;
    }


    // ═══════════════════════════════════════════════════════════════
    // GLMM LOGISTIQUE — PQL avec équations de Henderson
    // ═══════════════════════════════════════════════════════════════
    /**
     * Modèle logistique mixte à intercept aléatoire :
     *   logit(p_ij) = X_ij β + u_j ,  u_j ~ N(0, σ²_u)
     *
     * Algorithme : Penalized Quasi-Likelihood (PQL) de Breslow & Clayton
     * (1993) avec équations mixtes de Henderson :
     *   [ X'WX    X'WZ                 ] [β]   [X'Wz]
     *   [ Z'WX    Z'WZ + (1/σ²)·I_G    ] [u] = [Z'Wz]
     *
     * où z est la pseudo-réponse adjustée z_ij = η_ij + (y_ij−p_ij)/w_ij
     * et w_ij = p_ij(1−p_ij). σ² est mis à jour à chaque itération.
     *
     * ICC (intraclass correlation) = σ²_u / (σ²_u + π²/3) pour la lien logit.
     *
     * @param array $y         outcome binaire (n)
     * @param array $X         design matrix sans intercept (n×k)
     * @param array $cluster   identifiant cluster par observation (n)
     * @param array $names     noms des prédicteurs fixes
     * @param int   $max_iter  itérations max PQL (default 50)
     * @return array           résultats complets
     */
    public function glmmLogistic(array $y, array $X, array $cluster, array $names = [], int $max_iter = 50): array {
        $n = count($y);
        $k = count($X[0] ?? []);
        if($n < 30 || $k < 1) return ['error' => 'insufficient data'];

        // Mapping cluster_id → index 0..G-1
        $clu_map = [];
        $clu_list = array_values(array_unique($cluster));
        sort($clu_list);
        foreach($clu_list as $idx => $cid) $clu_map[$cid] = $idx;
        $G = count($clu_list);
        if($G < 3) return ['error' => 'need ≥ 3 clusters', 'n_clusters' => $G];
        $g_idx = array_map(fn($c) => $clu_map[$c], $cluster);

        // Construire matrice de design fixe (avec intercept)
        $Xf = [];
        for($i = 0; $i < $n; $i++){
            $row = [1.0];
            for($j = 0; $j < $k; $j++) $row[] = (float)$X[$i][$j];
            $Xf[] = $row;
        }
        $p = $k + 1; // nombre d'effets fixes (intercept + k)

        // Initialiser β par régression logistique standard
        $init = $this->logisticRegressionMulti($y, $X, $names);
        $beta = $init['converged'] ? $init['coef'] : array_fill(0, $p, 0.0);
        $u = array_fill(0, $G, 0.0);
        $sigma2 = 1.0; // valeur initiale
        $converged_pql = false;
        $iter_pql = 0;

        for($it = 0; $it < $max_iter; $it++){
            $iter_pql = $it + 1;
            // η_i = X_i β + u_{g(i)}, p_i, w_i, z_i
            $eta = array_fill(0, $n, 0.0);
            $w   = array_fill(0, $n, 0.0);
            $z   = array_fill(0, $n, 0.0);
            for($i = 0; $i < $n; $i++){
                $e = $u[$g_idx[$i]];
                for($j = 0; $j < $p; $j++) $e += $beta[$j] * $Xf[$i][$j];
                $e = max(min($e, 30), -30);
                $pi = 1.0 / (1.0 + exp(-$e));
                $wi = max(1e-6, $pi * (1 - $pi));
                $eta[$i] = $e;
                $w[$i]   = $wi;
                $z[$i]   = $e + ($y[$i] - $pi) / $wi;
            }

            // Construire les blocs d'Henderson (système (p+G) × (p+G))
            $M_size = $p + $G;
            $M = [];
            for($a = 0; $a < $M_size; $a++) $M[$a] = array_fill(0, $M_size, 0.0);
            $rhs = array_fill(0, $M_size, 0.0);

            // X'WX (p×p), X'WZ (p×G), X'Wz (p), Z'WZ diagonale (G×G), Z'Wz (G)
            for($i = 0; $i < $n; $i++){
                $wi = $w[$i]; $zi = $z[$i]; $gi = $g_idx[$i];
                // X'WX et X'Wz
                for($a = 0; $a < $p; $a++){
                    $xa_w = $Xf[$i][$a] * $wi;
                    for($b = 0; $b < $p; $b++) $M[$a][$b] += $xa_w * $Xf[$i][$b];
                    // X'WZ (bloc colonne droite-haut)
                    $M[$a][$p + $gi] += $xa_w;
                    // X'Wz
                    $rhs[$a] += $xa_w * $zi / $wi * $wi; // = X'·W·z
                }
                // Z'WX (bloc bas-gauche, symétrique de X'WZ)
                for($b = 0; $b < $p; $b++) $M[$p + $gi][$b] += $Xf[$i][$b] * $wi;
                // Z'WZ (diagonale)
                $M[$p + $gi][$p + $gi] += $wi;
                // Z'Wz
                $rhs[$p + $gi] += $wi * $zi;
            }
            // X'Wz : recompute corrigé
            for($a = 0; $a < $p; $a++){
                $s = 0.0;
                for($i = 0; $i < $n; $i++) $s += $Xf[$i][$a] * $w[$i] * $z[$i];
                $rhs[$a] = $s;
            }

            // Ajouter (1/σ²) à la diagonale du bloc bas-droit Z'WZ
            $inv_s2 = 1.0 / max(1e-6, $sigma2);
            for($j = 0; $j < $G; $j++) $M[$p + $j][$p + $j] += $inv_s2;

            // Résoudre M · [β; u] = rhs par inversion (taille p+G ≈ 26 pour 14 écoles+12 pred)
            $inv = $this->matrixInverse($M);
            if($inv === null) break;

            $sol = array_fill(0, $M_size, 0.0);
            for($a = 0; $a < $M_size; $a++){
                $s = 0.0;
                for($b = 0; $b < $M_size; $b++) $s += $inv[$a][$b] * $rhs[$b];
                $sol[$a] = $s;
            }

            $beta_new = array_slice($sol, 0, $p);
            $u_new    = array_slice($sol, $p, $G);

            // Mettre à jour σ² (estimateur ML : moyenne des u_j² + correction trace)
            // avec facteur de relaxation pour stabilité numérique
            $sum_u2 = 0.0;
            for($j = 0; $j < $G; $j++) $sum_u2 += $u_new[$j] ** 2;
            $trace_inv = 0.0;
            for($j = 0; $j < $G; $j++) $trace_inv += $inv[$p+$j][$p+$j];
            $sigma2_target = max(1e-4, ($sum_u2 + $trace_inv) / $G);
            // Relaxation 50% pour éviter les oscillations
            $sigma2_new = 0.5 * $sigma2 + 0.5 * $sigma2_target;

            // Convergence
            $maxd = 0.0;
            for($a = 0; $a < $p; $a++) $maxd = max($maxd, abs($beta_new[$a] - $beta[$a]));
            for($j = 0; $j < $G; $j++) $maxd = max($maxd, abs($u_new[$j] - $u[$j]));
            // Convergence : critère relatif pour σ² (plus stable numériquement)
            $rel_var = abs($sigma2_new - $sigma2) / max(0.01, $sigma2);

            $beta = $beta_new; $u = $u_new; $sigma2 = $sigma2_new;
            if($maxd < 5e-4 && $rel_var < 0.01){ $converged_pql = true; break; }
        }

        // SE des effets fixes à partir du bloc inv[1:p, 1:p]
        $se = array_fill(0, $p, null);
        if(isset($inv) && $inv !== null){
            for($j = 0; $j < $p; $j++){
                $v = $inv[$j][$j];
                $se[$j] = $v > 0 ? sqrt($v) : null;
            }
        }

        // OR + IC95 + p
        $or = []; $cilow = []; $cihigh = []; $pval = [];
        for($j = 0; $j < $p; $j++){
            $or[$j] = exp($beta[$j]);
            if($se[$j] !== null){
                $cilow[$j]  = exp($beta[$j] - 1.96 * $se[$j]);
                $cihigh[$j] = exp($beta[$j] + 1.96 * $se[$j]);
                $z_stat = $beta[$j] / $se[$j];
                $pval[$j] = 2 * (1 - $this->normalCDF(abs($z_stat)));
            } else {
                $cilow[$j] = $cihigh[$j] = $pval[$j] = null;
            }
        }

        // ICC (logit link → variance résiduelle latente = π²/3)
        $icc = $sigma2 / ($sigma2 + (M_PI * M_PI) / 3.0);
        // DEFF approché : 1 + (n_bar − 1) · ICC
        $n_bar = $n / $G;
        $deff = 1.0 + ($n_bar - 1.0) * $icc;

        $var_names = ['(Intercept)'];
        for($j = 0; $j < $k; $j++) $var_names[] = $names[$j] ?? "X$j";

        return [
            'converged'  => $converged_pql,
            'iter'       => $iter_pql,
            'n'          => $n,
            'n_clusters' => $G,
            'mean_cluster_size' => round($n_bar, 1),
            'names'      => $var_names,
            'coef'       => array_map(fn($v) => round($v, 4), $beta),
            'se'         => array_map(fn($v) => $v === null ? null : round($v, 4), $se),
            'or'         => array_map(fn($v) => round($v, 3), $or),
            'ci_low'     => array_map(fn($v) => $v === null ? null : round($v, 3), $cilow),
            'ci_high'    => array_map(fn($v) => $v === null ? null : round($v, 3), $cihigh),
            'p'          => array_map(fn($v) => $v === null ? null : round($v, 4), $pval),
            'sigma2_u'   => round($sigma2, 4),
            'icc'        => round($icc, 4),
            'deff'       => round($deff, 3),
            'random_effects' => array_combine($clu_list, array_map(fn($v) => round($v, 4), $u)),
        ];
    }

    // ═══════════════════════════════════════════════════════════════
    // GEE — GENERALIZED ESTIMATING EQUATIONS (Liang & Zeger 1986)
    // ═══════════════════════════════════════════════════════════════
    /**
     * Régression logistique GEE avec structure de corrélation EXCHANGEABLE
     * et variance robuste « sandwich » (Liang & Zeger 1986).
     *
     * Équation estimante :  Σ_i D_i' V_i^{-1} (y_i - μ_i) = 0
     * où   D_i = ∂μ_i/∂β = A_i X_i
     *      V_i = A_i^{1/2} R_i(α) A_i^{1/2}  avec  R_i(α) = (1-α)I + α11'
     *
     * Variance robuste : V_β = M_0^{-1} · M_1 · M_0^{-1}
     *   M_0 = Σ D_i' V_i^{-1} D_i  (information modèle)
     *   M_1 = Σ D_i' V_i^{-1} (y_i-μ_i)(y_i-μ_i)' V_i^{-1} D_i
     *
     * @param array $y         outcome binaire (n)
     * @param array $X         design matrix sans intercept (n×k)
     * @param array $cluster   identifiant cluster par observation (n)
     * @param array $names     noms des prédicteurs
     * @param int   $max_iter  itérations max (default 30)
     * @return array
     */
    public function geeLogistic(array $y, array $X, array $cluster, array $names = [], int $max_iter = 30): array {
        $n = count($y);
        $k = count($X[0] ?? []);
        if($n < 30 || $k < 1) return ['error' => 'insufficient data'];

        // Indexer les clusters
        $clu_map = [];
        $clu_list = array_values(array_unique($cluster));
        sort($clu_list);
        foreach($clu_list as $idx => $cid) $clu_map[$cid] = $idx;
        $G = count($clu_list);
        if($G < 3) return ['error' => 'need ≥ 3 clusters'];

        // Regrouper indices par cluster
        $groups = array_fill(0, $G, []);
        for($i = 0; $i < $n; $i++) $groups[$clu_map[$cluster[$i]]][] = $i;

        // Design avec intercept
        $Xf = [];
        for($i = 0; $i < $n; $i++){
            $row = [1.0];
            for($j = 0; $j < $k; $j++) $row[] = (float)$X[$i][$j];
            $Xf[] = $row;
        }
        $p = $k + 1;

        // Initialiser β par GLM logistique indépendant
        $init = $this->logisticRegressionMulti($y, $X, $names);
        $beta = $init['converged'] ? $init['coef'] : array_fill(0, $p, 0.0);
        $alpha = 0.0;
        $converged = false;
        $iter_gee = 0;

        for($it = 0; $it < $max_iter; $it++){
            $iter_gee = $it + 1;
            // 1. Calcul μ et résidus de Pearson
            $mu = array_fill(0, $n, 0.0);
            $a_diag = array_fill(0, $n, 0.0);
            $pr = array_fill(0, $n, 0.0);
            for($i = 0; $i < $n; $i++){
                $eta = 0.0;
                for($j = 0; $j < $p; $j++) $eta += $beta[$j] * $Xf[$i][$j];
                $eta = max(min($eta, 30), -30);
                $mui = 1.0 / (1.0 + exp(-$eta));
                $a_i = max(1e-6, $mui * (1 - $mui));
                $mu[$i] = $mui;
                $a_diag[$i] = $a_i;
                $pr[$i] = ($y[$i] - $mui) / sqrt($a_i);
            }

            // 2. Estimer α (exchangeable) : Σ_i Σ_{j<k} r_ij r_ik / [Σ_i n_i(n_i-1)/2 - p]
            $num = 0.0; $den_pairs = 0.0;
            foreach($groups as $g_indices){
                $ng = count($g_indices);
                if($ng < 2) continue;
                for($a = 0; $a < $ng; $a++) for($b = $a + 1; $b < $ng; $b++)
                    $num += $pr[$g_indices[$a]] * $pr[$g_indices[$b]];
                $den_pairs += $ng * ($ng - 1) / 2;
            }
            $alpha_new = ($den_pairs > $p) ? $num / ($den_pairs - $p) : 0.0;
            $alpha_new = max(-0.99, min(0.99, $alpha_new));

            // 3. Construire M0 = Σ D_i' V_i^{-1} D_i et S = Σ D_i' V_i^{-1} (y_i-μ_i)
            $M0 = []; for($a=0;$a<$p;$a++) $M0[$a]=array_fill(0,$p,0.0);
            $S  = array_fill(0, $p, 0.0);
            $M1 = []; for($a=0;$a<$p;$a++) $M1[$a]=array_fill(0,$p,0.0);

            foreach($groups as $g_indices){
                $ng = count($g_indices);
                if($ng === 0) continue;
                // X_i (ng×p), y_i, μ_i, A_i^{1/2}
                $X_i = []; $y_i = []; $mu_i = []; $a_sq = [];
                foreach($g_indices as $i){
                    $X_i[] = $Xf[$i]; $y_i[] = $y[$i]; $mu_i[] = $mu[$i];
                    $a_sq[] = sqrt($a_diag[$i]);
                }
                // R_i = (1-α) I + α 11', V_i = A^{1/2} R A^{1/2}
                // V_i^{-1} = A^{-1/2} R^{-1} A^{-1/2}
                // R^{-1} pour exchangeable : R = (1-α) I + α J ; J = 11'
                // R^{-1} = 1/(1-α) [I - α/(1-α+ng·α) · J]
                $denomR = (1 - $alpha_new) * (1 - $alpha_new + $ng * $alpha_new);
                $coef_I = 1.0 / (1 - $alpha_new);
                $coef_J = $denomR > 0 ? -$alpha_new / $denomR : 0.0;
                if($ng === 1){ $coef_I = 1.0; $coef_J = 0.0; }
                // V_i^{-1}[a][b] = (1/(a_sq[a]·a_sq[b])) · (coef_I·δ_ab + coef_J)
                $Vinv = [];
                for($a = 0; $a < $ng; $a++){
                    $row = array_fill(0, $ng, 0.0);
                    for($b = 0; $b < $ng; $b++){
                        $r_ab = ($a === $b) ? $coef_I + $coef_J : $coef_J;
                        $row[$b] = $r_ab / ($a_sq[$a] * $a_sq[$b]);
                    }
                    $Vinv[] = $row;
                }
                // D_i = A_i · X_i  (ng×p)
                $D_i = [];
                for($a = 0; $a < $ng; $a++){
                    $row = [];
                    foreach($X_i[$a] as $val) $row[] = $a_diag[$g_indices[$a]] * $val;
                    $D_i[] = $row;
                }
                // résidu r_i = y_i - μ_i (ng)
                $r_i = [];
                for($a = 0; $a < $ng; $a++) $r_i[] = $y_i[$a] - $mu_i[$a];

                // VinvD = V^-1 · D_i  (ng×p)
                $VinvD = [];
                for($a = 0; $a < $ng; $a++){
                    $row = array_fill(0, $p, 0.0);
                    for($b = 0; $b < $ng; $b++) for($j = 0; $j < $p; $j++)
                        $row[$j] += $Vinv[$a][$b] * $D_i[$b][$j];
                    $VinvD[] = $row;
                }
                // M0 += D_i' · VinvD
                for($a = 0; $a < $p; $a++) for($b = 0; $b < $p; $b++)
                    for($c = 0; $c < $ng; $c++) $M0[$a][$b] += $D_i[$c][$a] * $VinvD[$c][$b];
                // Vinv_r = V^-1 · r_i  (ng)
                $Vinv_r = array_fill(0, $ng, 0.0);
                for($a = 0; $a < $ng; $a++) for($b = 0; $b < $ng; $b++) $Vinv_r[$a] += $Vinv[$a][$b] * $r_i[$b];
                // S += D_i' · Vinv_r
                for($a = 0; $a < $p; $a++) for($c = 0; $c < $ng; $c++) $S[$a] += $D_i[$c][$a] * $Vinv_r[$c];
                // M1 += (D_i' Vinv r_i)(D_i' Vinv r_i)' (terme sandwich par cluster)
                $tmp = array_fill(0, $p, 0.0);
                for($a = 0; $a < $p; $a++) for($c = 0; $c < $ng; $c++) $tmp[$a] += $D_i[$c][$a] * $Vinv_r[$c];
                for($a = 0; $a < $p; $a++) for($b = 0; $b < $p; $b++) $M1[$a][$b] += $tmp[$a] * $tmp[$b];
            }

            // 4. Mettre à jour β : β_new = β + M0^{-1} S
            $invM0 = $this->matrixInverse($M0);
            if($invM0 === null) break;
            $delta = array_fill(0, $p, 0.0);
            for($a = 0; $a < $p; $a++) for($b = 0; $b < $p; $b++) $delta[$a] += $invM0[$a][$b] * $S[$b];

            $maxd = 0.0;
            for($a = 0; $a < $p; $a++){ $beta[$a] += $delta[$a]; $maxd = max($maxd, abs($delta[$a])); }
            $alpha = $alpha_new;
            if($maxd < 1e-5){ $converged = true; break; }
        }

        // Sandwich variance : V_β = invM0 · M1 · invM0
        $se_model = array_fill(0, $p, null);
        $se_robust = array_fill(0, $p, null);
        if(isset($invM0) && $invM0 !== null){
            for($j = 0; $j < $p; $j++){
                $vm = $invM0[$j][$j];
                $se_model[$j] = $vm > 0 ? sqrt($vm) : null;
            }
            // V_robust = invM0 · M1 · invM0
            $tmp = $this->matMul($invM0, $M1);
            $V_rob = $this->matMul($tmp, $invM0);
            for($j = 0; $j < $p; $j++){
                $vr = $V_rob[$j][$j];
                $se_robust[$j] = $vr > 0 ? sqrt($vr) : null;
            }
        }

        // OR + IC95 + p (utilise se_robust = sandwich)
        $or = []; $cilow = []; $cihigh = []; $pval = [];
        for($j = 0; $j < $p; $j++){
            $or[$j] = exp($beta[$j]);
            if($se_robust[$j] !== null){
                $cilow[$j]  = exp($beta[$j] - 1.96 * $se_robust[$j]);
                $cihigh[$j] = exp($beta[$j] + 1.96 * $se_robust[$j]);
                $z_stat = $beta[$j] / $se_robust[$j];
                $pval[$j] = 2 * (1 - $this->normalCDF(abs($z_stat)));
            } else {
                $cilow[$j] = $cihigh[$j] = $pval[$j] = null;
            }
        }

        $var_names = ['(Intercept)'];
        for($j = 0; $j < $k; $j++) $var_names[] = $names[$j] ?? "X$j";

        return [
            'converged'  => $converged,
            'iter'       => $iter_gee,
            'n'          => $n,
            'n_clusters' => $G,
            'correlation_structure' => 'exchangeable',
            'alpha'      => round($alpha, 4),
            'names'      => $var_names,
            'coef'       => array_map(fn($v) => round($v, 4), $beta),
            'se_model'   => array_map(fn($v) => $v === null ? null : round($v, 4), $se_model),
            'se_robust'  => array_map(fn($v) => $v === null ? null : round($v, 4), $se_robust),
            'or'         => array_map(fn($v) => round($v, 3), $or),
            'ci_low'     => array_map(fn($v) => $v === null ? null : round($v, 3), $cilow),
            'ci_high'    => array_map(fn($v) => $v === null ? null : round($v, 3), $cihigh),
            'p'          => array_map(fn($v) => $v === null ? null : round($v, 4), $pval),
        ];
    }


    // ═══════════════════════════════════════════════════════════════
    // MICE — MULTIPLE IMPUTATION BY CHAINED EQUATIONS (PMM)
    // ═══════════════════════════════════════════════════════════════
    /**
     * Imputation multiple par équations chaînées (van Buuren &
     * Groothuis-Oudshoorn 2011) avec Predictive Mean Matching (PMM).
     *
     * Algorithme :
     *   1. Init : remplir chaque valeur manquante par moyenne / mode
     *   2. Pour t = 1..max_iter (typiquement 20) :
     *      Pour chaque variable X_j avec des manquants :
     *        a. Sur les observations complètes : β_j = OLS/logit(X_j ~ X_{-j})
     *        b. Pour chaque cas manquant i :
     *             ŷ_i = X_{-j,i} β_j
     *           Trouver les d=5 cas observés dont ŷ est le plus proche de ŷ_i
     *           Tirer aléatoirement parmi ces d donneurs et copier sa valeur
     *   3. Répéter m fois pour produire m datasets imputés
     *
     * @param array $data       n × p tableau (NULL = manquant)
     * @param array $var_types  types : 'continuous' | 'binary' | 'categorical'
     * @param int   $m          nombre d'imputations
     * @param int   $max_iter   itérations par imputation (default 20)
     * @param int   $donors     candidats PMM (default 5)
     * @return array            ['imputations' => [...$m datasets...], 'meta' => ...]
     */
    public function mice(array $data, array $var_types, int $m = 20, int $max_iter = 20, int $donors = 5): array {
        $n = count($data);
        if($n === 0) return ['error' => 'empty data'];
        $p = count($data[0]);
        $var_names = array_keys($data[0]);

        // 1. Détecter le pattern de manquants
        $is_missing = [];
        $n_missing_per_var = array_fill_keys($var_names, 0);
        for($i = 0; $i < $n; $i++){
            $row = [];
            foreach($var_names as $v){
                $val = $data[$i][$v] ?? null;
                $miss = ($val === null || $val === '' || (is_numeric($val) && is_nan((float)$val)));
                $row[$v] = $miss;
                if($miss) $n_missing_per_var[$v]++;
            }
            $is_missing[] = $row;
        }
        $vars_to_impute = array_filter($var_names, fn($v) => $n_missing_per_var[$v] > 0);

        // 2. Imputation initiale (moyenne / mode)
        $initial_values = [];
        foreach($var_names as $v){
            $vals = [];
            for($i = 0; $i < $n; $i++) if(!$is_missing[$i][$v]) $vals[] = $data[$i][$v];
            if(count($vals) === 0){ $initial_values[$v] = 0; continue; }
            if(($var_types[$v] ?? 'continuous') === 'continuous'){
                $initial_values[$v] = array_sum($vals) / count($vals);
            } else {
                // mode
                $counts = array_count_values(array_map('strval', $vals));
                arsort($counts);
                $initial_values[$v] = array_key_first($counts);
            }
        }

        // 3. Construire m imputations
        $imputations = [];
        for($mm = 0; $mm < $m; $mm++){
            // Copier les données + init des manquants
            $imp = [];
            for($i = 0; $i < $n; $i++){
                $row = [];
                foreach($var_names as $v){
                    $row[$v] = $is_missing[$i][$v] ? $initial_values[$v] : $data[$i][$v];
                }
                $imp[] = $row;
            }

            // 4. Itérations Gibbs (chaînées)
            for($it = 0; $it < $max_iter; $it++){
                foreach($vars_to_impute as $j_var){
                    // Construire (y_obs, X_obs) et (X_mis) pour la variable j_var
                    $other_vars = array_diff($var_names, [$j_var]);
                    $obs_idx = []; $mis_idx = [];
                    for($i = 0; $i < $n; $i++) $is_missing[$i][$j_var] ? $mis_idx[] = $i : $obs_idx[] = $i;
                    if(count($obs_idx) < 10 || count($mis_idx) === 0) continue;

                    // Construire matrices
                    $X_obs = []; $y_obs = [];
                    foreach($obs_idx as $i){
                        $row = [1.0];
                        foreach($other_vars as $v) $row[] = (float)$imp[$i][$v];
                        $X_obs[] = $row;
                        $y_obs[] = (float)$imp[$i][$j_var];
                    }
                    $X_mis = [];
                    foreach($mis_idx as $i){
                        $row = [1.0];
                        foreach($other_vars as $v) $row[] = (float)$imp[$i][$v];
                        $X_mis[] = $row;
                    }

                    // Régression : OLS ou logistique selon le type
                    $type = $var_types[$j_var] ?? 'continuous';
                    if($type === 'binary'){
                        // logistique multivariée
                        $X_obs_noint = array_map(fn($r) => array_slice($r, 1), $X_obs);
                        $fit = $this->logisticRegressionMulti($y_obs, $X_obs_noint, []);
                        if(isset($fit['error']) || !$fit['converged']) continue;
                        $beta = $fit['coef'];
                    } else {
                        $reg = $this->olsRegression($y_obs, $X_obs);
                        if($reg === null) continue;
                        $beta = $reg['beta'];
                    }

                    // Prédictions sur observed + missing
                    $K = count($beta);
                    $yhat_obs = []; $yhat_mis = [];
                    foreach($X_obs as $r){ $e = 0; for($a=0;$a<$K;$a++) $e += $r[$a]*$beta[$a]; $yhat_obs[] = $e; }
                    foreach($X_mis as $r){ $e = 0; for($a=0;$a<$K;$a++) $e += $r[$a]*$beta[$a]; $yhat_mis[] = $e; }

                    // 5. PMM : pour chaque cas manquant, trouver d donneurs les plus proches
                    for($mi = 0; $mi < count($mis_idx); $mi++){
                        $target = $yhat_mis[$mi];
                        // Calculer distances absolues à tous les observed
                        $dists = [];
                        foreach($yhat_obs as $oi => $oh) $dists[$oi] = abs($oh - $target);
                        asort($dists);
                        $top = array_slice($dists, 0, $donors, true);
                        $donor_keys = array_keys($top);
                        $picked = $donor_keys[mt_rand(0, count($donor_keys) - 1)];
                        $value_to_use = $y_obs[$picked];
                        $imp[$mis_idx[$mi]][$j_var] = $value_to_use;
                    }
                }
            }
            $imputations[] = $imp;
        }

        return [
            'm'             => $m,
            'iterations'    => $max_iter,
            'donors_pmm'    => $donors,
            'n_total'       => $n,
            'vars_imputed'  => array_values($vars_to_impute),
            'pct_missing'   => array_map(fn($cnt) => round($cnt * 100 / $n, 2), $n_missing_per_var),
            'imputations'   => $imputations, // tableau de m jeux complets
        ];
    }

    /**
     * Combine les résultats de m analyses par les règles de Rubin (1987).
     *
     *   β̂  = mean(β̂_m)
     *   U  = mean(SE_m²)                  variance intra-imputation
     *   B  = (1/(m-1)) Σ (β̂_m - β̂)²       variance inter-imputation
     *   T  = U + (1 + 1/m) · B              variance totale
     *   SE = √T
     *   df = (m-1) · (1 + U / ((1+1/m)·B))²
     *
     * @param array $estimates    [[beta_1, beta_2, ...], ...] (m × p)
     * @param array $standard_err [[se_1, se_2, ...], ...]      (m × p)
     * @return array              ['beta'=>, 'se'=>, 'p'=>, 'ci_low'=>, 'ci_high'=>, 'df'=>, 'or'=>, ...]
     */
    public function rubinPool(array $estimates, array $standard_err): array {
        $m = count($estimates);
        if($m < 2) return ['error' => 'need at least 2 imputations'];
        $p = count($estimates[0]);

        $beta_bar = array_fill(0, $p, 0.0);
        for($i = 0; $i < $m; $i++) for($j = 0; $j < $p; $j++) $beta_bar[$j] += $estimates[$i][$j];
        for($j = 0; $j < $p; $j++) $beta_bar[$j] /= $m;

        // Within (U) et Between (B)
        $U = array_fill(0, $p, 0.0);
        $B = array_fill(0, $p, 0.0);
        for($i = 0; $i < $m; $i++) for($j = 0; $j < $p; $j++)
            $U[$j] += ($standard_err[$i][$j] ?? 0) ** 2;
        for($j = 0; $j < $p; $j++) $U[$j] /= $m;

        for($i = 0; $i < $m; $i++) for($j = 0; $j < $p; $j++)
            $B[$j] += ($estimates[$i][$j] - $beta_bar[$j]) ** 2;
        for($j = 0; $j < $p; $j++) $B[$j] /= ($m - 1);

        // Variance totale T et SE
        $T = []; $SE = []; $df = []; $p_val = []; $cilo = []; $cihi = []; $or = []; $or_lo = []; $or_hi = [];
        for($j = 0; $j < $p; $j++){
            $T_j = $U[$j] + (1 + 1.0/$m) * $B[$j];
            $T[$j]   = $T_j;
            $SE[$j]  = $T_j > 0 ? sqrt($T_j) : null;
            // df Barnard & Rubin (1999) simplifié
            $r_j = (1 + 1.0/$m) * $B[$j] / max(1e-12, $U[$j]);
            $df_j = ($m - 1) * (1 + 1.0/$r_j) ** 2;
            $df[$j] = $df_j;

            if($SE[$j] !== null){
                $z = $beta_bar[$j] / $SE[$j];
                $p_val[$j] = 2 * (1 - $this->normalCDF(abs($z)));
                $cilo[$j] = $beta_bar[$j] - 1.96 * $SE[$j];
                $cihi[$j] = $beta_bar[$j] + 1.96 * $SE[$j];
                $or[$j]    = exp($beta_bar[$j]);
                $or_lo[$j] = exp($cilo[$j]);
                $or_hi[$j] = exp($cihi[$j]);
            } else {
                $p_val[$j] = null; $cilo[$j] = null; $cihi[$j] = null;
                $or[$j] = null; $or_lo[$j] = null; $or_hi[$j] = null;
            }
        }

        return [
            'm'         => $m,
            'beta'      => array_map(fn($v) => round($v, 4), $beta_bar),
            'se'        => array_map(fn($v) => $v === null ? null : round($v, 4), $SE),
            'or'        => array_map(fn($v) => $v === null ? null : round($v, 3), $or),
            'ci_low'    => array_map(fn($v) => $v === null ? null : round($v, 3), $or_lo),
            'ci_high'   => array_map(fn($v) => $v === null ? null : round($v, 3), $or_hi),
            'p'         => array_map(fn($v) => $v === null ? null : round($v, 4), $p_val),
            'df'        => array_map(fn($v) => round($v, 1), $df),
            'within_U'  => array_map(fn($v) => round($v, 6), $U),
            'between_B' => array_map(fn($v) => round($v, 6), $B),
            'total_T'   => array_map(fn($v) => round($v, 6), $T),
            'fmi'       => array_map(function($U_j, $B_j) {  // fraction missing info
                $rT = $U_j + (1 + 1.0/count(func_get_args())) * $B_j;
                return round((1 + 1.0/count(func_get_args())) * $B_j / max(1e-12, $rT), 4);
            }, $U, $B),
        ];
    }


}
