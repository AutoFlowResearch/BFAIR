import React, { useContext } from 'react';
import { Study } from '../../../../../shared/model/investigation.model';
import InvestigationFormAssay from '../InvestigationFormAssay';
import addIcon from '../../../../../../assets/images/add-icon.svg';
import { addAssay, removeAssay } from '../../../../../shared/redux/actions';
import { connect } from 'react-redux';

const InvestigationFormAssays = (props) => {
  const { study, removeAssay, addAssay } = props;
  const { assays } = study;
  return (
    <React.Fragment>
      {assays.map((assay, index) => (
        <InvestigationFormAssay
          key={index}
          assay={assay}
          removeAssay={() => removeAssay(study, index)}
        />
      ))}
      <div className="section__add">
        <button
          className="section__add__button"
          onClick={() => addAssay(study)}
        >
          <img alt="" src={addIcon} className="section__add__button__icon" />
          Add New Assay
        </button>
      </div>
    </React.Fragment>
  );
};

const mapStateToProps = (state) => ({});

const mapDispatchToProps = (dispatch) => ({
  removeAssay: (study, index) => dispatch(removeAssay(study, index)),
  addAssay: (study) => dispatch(addAssay(study)),
});

export default connect(
  mapStateToProps,
  mapDispatchToProps,
)(InvestigationFormAssays);
