# Invio a Statistics in Medicine: cosa c'è e cosa resta da fare

Pacchetto preparato il 18 novembre 2026 dalla versione 9 inviata a SMMR (zip del 10 novembre). Il testo scientifico non è cambiato: cambiano solo le parti che dipendono dalla rivista.

## Cosa caricare e dove (sistema Research Exchange, submission.wiley.com)

| File | Tipo di file nel sistema |
|---|---|
| 01_cover_letter/MMM_SiM_cover_letter.docx | Cover letter |
| 02_manuscript/MMM_SiM_manuscript.docx | Main document |
| MMM_SiM_review_copy.pdf | PDF unico per la revisione (manoscritto e supplemento, 47 pagine): la rivista lo consiglia per il formato libero |
| 03_figures/Figure1-3 (.tif a 300 dpi, oppure .pdf) | Figure |
| 04_supporting_information/MMM_SiM_supporting_information.docx (o .pdf) | Supporting Information |
| 05_data_files/fewclust_0.4.2_source.tar.gz | Data Files (sorgente del pacchetto fewclust 0.4.2, tag 0.4.2 del repository pubblico) |

Il repository del codice (434 MB) e i risultati per replicato (9,6 GB) sono troppo grandi da caricare: restano su Zenodo e sono citati come riferimenti 35-37.

## Cosa è cambiato rispetto alla versione per SMMR

1. Revisione non anonima per gli autori (single-anonymous): autori, affiliazioni, autore corrispondente, running head e conteggio delle parole sono ora nella prima pagina del manoscritto. Non serve una title page separata.
2. Dichiarazioni in fondo al manoscritto, prima della bibliografia: Data availability statement, Author contributions, Funding, Conflict of interest, Ethics statement. Il testo è quello della title page v9.
3. Dichiarazione sull'uso dell'intelligenza artificiale generativa spostata nei Metodi, in coda alla sezione 2.4 (Implementation), come chiede la rivista; il testo non cambia. Il conteggio del testo principale sale di 83 parole, da 4.676 a 4.759.
4. Dichiarazione sui dati con i link espliciti a GitHub e i tre depositi Zenodo (codice release 1.1, fewclust 0.4.2, risultati per replicato) citati come riferimenti 35-37, perché la rivista chiede di citare i dati come gli articoli.
5. Parole chiave da 7 a 6 (massimo della rivista): tolta "simulation study".
6. Supplemento intitolato "Supporting Information", con titolo e autori in testa. Larghezze delle colonne fissate perché le tabelle non si spezzino nel PDF.
7. Tabella 5 del manoscritto: larghezze delle colonne fissate; prima la conversione in PDF spezzava le parole.
8. Lettera all'editor nuova, indirizzata agli editor di Statistics in Medicine.

Restano invariati: titolo, abstract (246 parole, limite 250), il resto del testo da Introduzione a Conclusioni, Tabelle 1-4, Figure 1-3, riferimenti 1-34.

## Da fare prima di premere Submit

1. **Release del codice: fatto.** Release 1.1 di clustered_paper (tag 1.1, DOI 10.5281/zenodo.24345879), con il codice della riscrittura v9 (py/v9, R/v9, tables/v9, manuscript/v9; rieseguito il 18 novembre, riproduce identici i file della v9). Riferimento 35 e dichiarazione sui dati la citano.
2. ORCID di tutti gli autori ed email dei coautori: il sistema li chiede.
3. Data della lettera (ora 18 novembre 2026).
4. Se il questionario chiede di invii precedenti, il lavoro è stato inviato due volte a Statistical Methods in Medical Research (SMM-26-1128).
5. Revisori suggeriti: facoltativi, non indicati.
6. Rileggere la dichiarazione sull'uso dell'intelligenza artificiale (sezione 2.4) e confermare che è ancora esatta.
7. Il riassunto non tecnico (lay abstract) si manda solo dopo l'accettazione.
