import { Assay, Study } from '../model/investigation.model';
import {
  ADD_ASSAY,
  ADD_STUDY,
  REMOVE_ASSAY,
  REMOVE_STUDY,
  SET_ASSAY_TITLE,
  SET_STUDY_TITLE,
} from '../../components/templates/InvestigationFormPage/InvestigationFormPage.constants';

export const addStudy = () => ({
  type: ADD_STUDY,
});

export const removeStudy = (index: number) => ({
  type: REMOVE_STUDY,
  payload: index,
});

export const addAssay = (study: Study) => ({
  type: ADD_ASSAY,
  payload: { study },
});

export const removeAssay = (study: Study, index: number) => ({
  type: REMOVE_ASSAY,
  payload: { study, index },
});

export const setStudyTitle = (study: Study, title: string) => ({
  type: SET_STUDY_TITLE,
  payload: { study, title },
});

export const setAssayTitle = (assay: Assay, title: string) => ({
  type: SET_ASSAY_TITLE,
  payload: { assay, title },
});
