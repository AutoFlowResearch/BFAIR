import React from 'react';
import { Assay } from '../../../../../shared/model/investigation.model';
import Collapsible from '../../../../atoms/Collapsible';
import DatePickerInput from '../../../../atoms/DatePickerInput';
import DeleteIcon from '../../../../atoms/DeleteIcon';
import Input from '../../../../atoms/Input';

import './InvestigationFormAssay.style.scss';

import fileIcon from '../../../../../../assets/images/file-icon.svg';
import addIcon from '../../../../../../assets/images/add-icon.svg';
import { setAssayTitle } from '../../../../../shared/redux/actions';
import { connect } from 'react-redux';
import {
  MEASUREMENT_TYPE_ANNOTATIONS,
  PERFORMERS,
  TECHNOLOGY_TYPE_ANNOTATIONS,
} from '../../../../../shared/constants';
import SelectBar from '../../../../atoms/SelectBar';

const InvestigationFormAssayTitle = (props) => {
  const { assay, removeAssay } = props;
  return (
    <div className="assay__title">
      {assay.title || 'New Assay'}
      <div className="assay__delete">
        <DeleteIcon onClick={removeAssay} />
      </div>
    </div>
  );
};

const InvestigationFormAssayFile = (props) => {
  return (
    <div className="assay__file">
      <div className="row">
        <div className="col-sm-12 col-lg-6 col-xl-4">
          <p className="assay__file__object-title">Digital Object</p>
          <div className="assay__file__upload-btn">
            <img alt="" src={fileIcon} />
            File Upload
          </div>
        </div>
        <div className="col-sm-12 col-lg-6 col-xl-4">
          <Input label="Bio Material Object" />
        </div>
      </div>
      <div className="section__add">
        <button className="section__add__button">
          <img alt="" src={addIcon} className="section__add__button__icon" />
          Add New File
        </button>
      </div>
    </div>
  );
};

const InvestigationFormAssay = (props: {
  assay: Assay;
  removeAssay: any;
  setAssayTitle?: any;
}) => {
  const { assay, removeAssay } = props;

  const handleAssayTitleChange = (title) => {
    setAssayTitle(assay, title);
  };
  return (
    <div className="investigation-form__assay">
      <Collapsible
        trigger={
          <InvestigationFormAssayTitle
            assay={assay}
            removeAssay={removeAssay}
          />
        }
      >
        <Input label="Title" onChange={handleAssayTitleChange} />
        <div className="row">
          <div className="col-sm-12 col-lg-6 col-xl-4">
            <DatePickerInput label="Start Date" />
          </div>
          <div className="col-sm-12 col-lg-6 col-xl-4">
            <DatePickerInput label="End Date" />
          </div>
          <div className="col-sm-12 col-lg-6 col-xl-4">
            <Input label="Run Order" />
          </div>
        </div>
        <SelectBar
          label={'Performer'}
          placeHolderText={'Performer'}
          options={PERFORMERS}
        />
        <Collapsible trigger="Inputs">
          <InvestigationFormAssayFile />
        </Collapsible>
        <Collapsible trigger="Outputs">
          <InvestigationFormAssayFile />
        </Collapsible>

        <h3 className="form-section__title">Measurement Type</h3>
        <Input label="Name" />
        <SelectBar
          label={'Annotation'}
          placeHolderText={'Annotation'}
          options={MEASUREMENT_TYPE_ANNOTATIONS}
        />
        <h3 className="form-section__title">Technology Type</h3>
        <Input label="Name" />
        <SelectBar
          label={'Annotation'}
          placeHolderText={'Annotation'}
          options={TECHNOLOGY_TYPE_ANNOTATIONS}
        />
      </Collapsible>
    </div>
  );
};

const mapDispatchToProps = (dispatch) => ({
  setAssayTitle: (assay, title) => dispatch(setAssayTitle(assay, title)),
});

export default connect(null, mapDispatchToProps)(InvestigationFormAssay);
